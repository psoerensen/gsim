#include "simulation_workflow.h"
#include "gaussian_draws.h"
#include "breeding_rng.h"
#include "meiosis_core.h"
#include "packed_chromosome.h"
#include "bed_storage.h"
#include "metadata_storage.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <locale>
#include <map>
#include <set>
#include <numeric>
#include <algorithm>
#include <limits>
namespace gsim::native {
namespace {
using U=std::uint64_t;using V=std::vector<double>;namespace fs=std::filesystem;
void need(bool x,const std::string& s){if(!x)throw std::runtime_error(s);}
std::string id(const std::string& s){need(!s.empty()&&s!="0"&&s.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-")==s.npos,"workflow identities must use nonempty ASCII letters/digits/_.- and cannot be 0");return s;}
std::vector<std::string> words(const std::string&s){std::istringstream in(s);in.imbue(std::locale::classic());std::vector<std::string> out;std::string x;while(in>>x)out.push_back(x);return out;}
U count(const std::string&s){need(!s.empty()&&s.find_first_not_of("0123456789")==s.npos,"integer setting must contain decimal digits");std::size_t at=0;auto n=std::stoull(s,&at);need(at==s.size(),"invalid integer");return n;}
double real(const std::string&s){std::istringstream in(s);in.imbue(std::locale::classic());double x;in>>x;need(bool(in)&&std::isfinite(x),"nonfinite/invalid numeric setting");in>>std::ws;need(in.eof(),"invalid numeric suffix");return x;}
V factor(const V&a,U q){return gaussian_lower(a,q);}

std::ofstream output(const fs::path&p){std::ofstream f(p,std::ios::binary);f.imbue(std::locale::classic());f<<std::setprecision(17);need(bool(f),"cannot open output: "+p.string());return f;}
U checksum(const fs::path&p){std::ifstream f(p,std::ios::binary);need(bool(f),"cannot checksum input");U h=14695981039346656037ULL;char buf[65536];while(f){f.read(buf,sizeof buf);for(std::streamsize i=0;i<f.gcount();++i){h^=static_cast<unsigned char>(buf[i]);h*=1099511628211ULL;}}need(f.eof(),"checksum read error");return h;}
std::vector<U> sample(const std::vector<metadata::SimulationSampleMetadata>&p,U k,U seed){std::vector<U> ix(p.size()),keys(p.size());std::iota(ix.begin(),ix.end(),0);for(U i=0;i<p.size();++i)keys[i]=breeding::stream(seed,p[i].individual_id,"rank","sampling")();std::partial_sort(ix.begin(),ix.begin()+k,ix.end(),[&](U a,U b){return keys[a]==keys[b]?p[a].individual_id<p[b].individual_id:keys[a]<keys[b];});ix.resize(k);return ix;}
}
void simulate_recipe(const std::string& wire,const std::string& destination){
 if(wire.rfind("GSIM_GAUSSIAN 1",0)==0){simulate_gaussian_recipe(wire,destination);return;}
 std::istringstream config(wire);config.imbue(std::locale::classic());std::string magic;std::getline(config,magic);if(!magic.empty()&&magic.back()=='\r')magic.pop_back();need(magic=="GSIM_SIMULATION 1","unsupported simulation wire version");
 std::map<std::string,std::string> c;std::string key,value;while(config>>key){need(bool(config>>std::quoted(value)),"malformed simulation setting");need(c.emplace(key,value).second,"duplicate simulation setting");}need(config.eof(),"malformed simulation settings");
 const std::set<std::string> allowed{"pedigree","variants","founders","effects","traits","seed","genotyped","phenotyped","records_per_individual","visits","missing_rate","residual_covariance","maximum_memory_bytes","export_truth_phases"};for(const auto&kv:c)need(allowed.count(kv.first),"unknown simulation setting: "+kv.first);
 for(auto k:{"pedigree","variants","founders","effects","traits"})need(c.count(k),std::string("missing simulation setting: ")+k);
 auto get=[&](const char*k,const char*d){return c.count(k)?c.at(k):std::string(d);};
 auto traits=words(c.at("traits"));U q=traits.size();need(q>0&&q<=64,"workflow supports 1..64 traits");std::set<std::string> distinct;for(auto&t:traits)need(distinct.insert(id(t)).second,"duplicate trait");
 need(!(c.count("records_per_individual")&&c.count("visits")),"specify records_per_individual or legacy visits, not both");
 U seed=count(get("seed","1")),visits=count(c.count("records_per_individual")?c.at("records_per_individual"):get("visits","1")),budget=count(get("maximum_memory_bytes","4294967296"));need(seed<=9007199254740981ULL&&visits>0,"invalid workflow seed or records_per_individual");double missing=real(get("missing_rate","0"));need(missing>=0&&missing<1,"missing_rate must be in [0,1)");
 std::map<std::string,U> input_identity;for(auto k:{"pedigree","variants","founders","effects"})input_identity[k]=checksum(c.at(k));
 auto pedigree=metadata::ValidatedSampleMetadata::read_fam(c.at("pedigree"));const auto& pop=pedigree.records();U n=pop.size();auto variants=metadata::ValidatedVariantMetadata::read_bim(c.at("variants"));U m=variants.size();need(n>0&&m>0,"empty simulation domain");
 U ng=c.count("genotyped")?count(c.at("genotyped")):n,np=c.count("phenotyped")?count(c.at("phenotyped")):n;need(ng>0&&ng<=n&&np>0&&np<=n,"selection count out of range");
 long double bytes=16.L*((n+63)/64)*m+16.L*q*n+16.L*q*m+16.L*ng*std::min<U>(m,128)+256.L*n;
 need(bytes<=budget,"simulation working-memory admission exceeds maximum_memory_bytes");need(n<=std::numeric_limits<std::size_t>::max()/q&&m<=std::numeric_limits<std::size_t>::max()/q,"simulation extent overflow");
 V covariance(q*q,0);for(U t=0;t<q;++t)covariance[t*q+t]=1;
 if(c.count("residual_covariance")){auto w=words(c.at("residual_covariance"));need(w.size()==q*q,"residual covariance extent");for(U i=0;i<w.size();++i)covariance[i]=real(w[i]);}V L=factor(covariance,q);
 std::map<std::string,U> index;std::vector<U> sire(n,n),dam(n,n);for(U i=0;i<n;++i){id(pop[i].individual_id);need(index.emplace(pop[i].individual_id,i).second,"duplicate animal ID");}
 for(U i=0;i<n;++i){for(auto pair:{std::make_pair(pop[i].paternal_id,&sire[i]),std::make_pair(pop[i].maternal_id,&dam[i])})if(pair.first!="0"){need(index.count(pair.first)&&index.at(pair.first)<i,"workflow pedigree must be parent-first with both parents recorded");*pair.second=index.at(pair.first);}need((sire[i]==n)==(dam[i]==n),"workflow requires zero or two known parents");}
 V positions(m);std::set<std::string> marker_ids;auto chromosome=variants.records()[0].chromosome;
 for(U j=0;j<m;++j){const auto&v=variants.records()[j];id(v.marker_id);need(marker_ids.insert(v.marker_id).second,"duplicate marker");need(v.chromosome==chromosome,"initial workflow supports one chromosome");positions[j]=v.genetic_distance_cm*.01;need(j==0||positions[j]>=positions[j-1],"map must be nondecreasing");}
 PhasedHaplotypeMatrix h1(n,m),h2(n,m);std::set<U> supplied;std::ifstream founders(c.at("founders"));need(bool(founders),"cannot read founder phases");std::string animal,a,b;
 while(founders>>animal){need(bool(founders>>a>>b),"truncated founder phase row");need(index.count(animal),"unknown founder");U i=index.at(animal);need(sire[i]==n&&supplied.insert(i).second&&a.size()==m&&b.size()==m,"invalid founder phase row");for(U j=0;j<m;++j){need((a[j]=='0'||a[j]=='1')&&(b[j]=='0'||b[j]=='1'),"founder phases must be 0/1");h1.set_allele(i,j,a[j]-'0');h2.set_allele(i,j,b[j]-'0');}}need(founders.eof(),"truncated founder phases");for(U i=0;i<n;++i)if(sire[i]==n)need(supplied.count(i),"founder phases missing");
 for(U i=0;i<n;++i)if(sire[i]!=n){for(int side=1;side<=2;++side){auto e=meiosis::draw_events_raw(positions.data(),m,seed,breeding::hash(pop[i].individual_id),breeding::hash(chromosome),side);std::vector<U> boundaries;for(double at:e.crossovers)boundaries.push_back(std::lower_bound(positions.begin(),positions.end(),at)-positions.begin());(side==1?h1:h2).make_gamete(i,h1,h2,side==1?sire[i]:dam[i],e.starting_haplotype-1,boundaries.data(),boundaries.size());}}
 V effects(m*q);std::ifstream ef(c.at("effects"));need(bool(ef),"cannot read marker effects");for(U j=0;j<m;++j){std::string name;ef>>name;need(bool(ef)&&name==variants.records()[j].marker_id,"marker effect order mismatch");for(U t=0;t<q;++t){ef>>effects[j+t*m];need(bool(ef)&&std::isfinite(effects[j+t*m]),"nonfinite/absent marker effect");}}ef>>std::ws;need(ef.eof(),"extra marker effects");
 V centers(m,0),genetic(n*q,0);U founder_count=supplied.size();for(U j=0;j<m;++j){for(U i:supplied)centers[j]+=h1.allele(i,j)+h2.allele(i,j);centers[j]/=founder_count;for(U i=0;i<n;++i)for(U t=0;t<q;++t)genetic[i+t*n]+=(h1.allele(i,j)+h2.allele(i,j)-centers[j])*effects[j+t*m];}
 auto selected=sample(pop,ng,seed+1),observed=sample(pop,np,seed+2);std::sort(selected.begin(),selected.end());
 fs::path out=fs::absolute(destination).lexically_normal(),stage=out.string()+".staging";need(!fs::exists(out)&&!fs::exists(stage),"preserve existing simulation output/staging directory");need(fs::is_directory(out.parent_path()),"simulation output parent must exist");fs::create_directory(stage);
 if(c.count("export_truth_phases")) {
  need(c.at("export_truth_phases")=="0" || c.at("export_truth_phases")=="1","invalid export_truth_phases");
  if(c.at("export_truth_phases")=="1") {
   auto f=output(stage/"truth_phase_pairs.u8");
   for(U j=0;j<m;++j)for(U i=0;i<n;++i)f.put(char(h1.allele(i,j)+2*h2.allele(i,j)));
   f.close();need(bool(f),"truth phase export failed");
  }
 }
 // Existing packed BED writer; never materialize the dense genotype matrix.
 PhasedBedSink sink((stage/"genotypes.bed").string(),ng,false,128);
 for(U first=0;first<m;first+=128){U width=std::min<U>(128,m-first);PhasedHaplotypeMatrix a(ng,width),b(ng,width);for(U j=0;j<width;++j)for(U k=0;k<ng;++k){a.set_allele(k,j,h1.allele(selected[k],first+j));b.set_allele(k,j,h2.allele(selected[k],first+j));}sink.append(a,b);}sink.finalize();(void)variants.write_bim((stage/"genotypes.bim").string(),false);
 std::vector<metadata::SimulationSampleMetadata> chosen;for(U i:selected)chosen.push_back(pop[i]);(void)metadata::ValidatedSampleMetadata(chosen,true).write_fam((stage/"genotypes.fam").string(),false);
 auto ped=output(stage/"pedigree.csv");ped<<"animal,sire,dam\n";for(const auto&p:pop)ped<<p.individual_id<<','<<p.paternal_id<<','<<p.maternal_id<<'\n';ped.close();need(bool(ped),"pedigree write failed");
 auto truth=output(stage/"true_breeding_values.csv");truth<<"animal,trait,true_breeding_value\n";for(U i=0;i<n;++i)for(U t=0;t<q;++t)truth<<pop[i].individual_id<<','<<traits[t]<<','<<genetic[i+t*n]<<'\n';truth.close();need(bool(truth),"truth write failed");
 auto values=output(stage/"phenotypes.csv"),record=output(stage/"record_truth.csv");values<<"record_id,animal";for(const auto&t:traits)values<<','<<t;values<<'\n';record<<"record_id,animal,trait,genetic,residual,phenotype,observed\n";
 U blocks=0,scalars=0;for(U i:observed)for(U visit=0;visit<visits;++visit){auto unit=pop[i].individual_id+":"+std::to_string(visit+1);std::vector<bool> keep(q);V residual(q,0);for(U t=0;t<q;++t){keep[t]=breeding::stream(seed+4,unit,traits[t],"missingness").uniform()>=missing;for(U j=0;j<q;++j){auto e=breeding::stream(seed+3,unit,traits[j],"residual");residual[t]+=L[t+j*q]*e.normal();}}bool any=std::any_of(keep.begin(),keep.end(),[](bool x){return x;});if(any){values<<unit<<','<<pop[i].individual_id;++blocks;}
 for(U t=0;t<q;++t){double y=genetic[i+t*n]+residual[t];need(std::isfinite(y),"nonfinite phenotype");record<<unit<<','<<pop[i].individual_id<<','<<traits[t]<<','<<genetic[i+t*n]<<','<<residual[t]<<','<<y<<','<<(keep[t]?1:0)<<'\n';if(any){values<<',';if(keep[t]){values<<y;++scalars;}else values<<"NA";}}if(any)values<<'\n';}
 need(blocks>0,"all phenotypes missing");values.close();record.close();need(bool(values)&&bool(record),"record write failed");
 auto settings=output(stage/"effective_settings.json");settings<<"{\n  \"schema\": \"gsim.simulation/v1\",\n  \"seed\": "<<seed<<",\n  \"genotyped\": "<<ng<<",\n  \"phenotyped\": "<<np<<",\n  \"records_per_individual\": "<<visits<<",\n  \"visits\": "<<visits<<",\n  \"missing_rate\": "<<missing<<",\n  \"maximum_memory_bytes\": "<<budget<<",\n  \"traits\": [";for(U t=0;t<q;++t){if(t)settings<<',';settings<<std::quoted(traits[t]);}settings<<"],\n  \"residual_covariance_row_major\": [";for(U i=0;i<q*q;++i){if(i)settings<<',';settings<<covariance[i];}settings<<"]\n}\n";settings.close();need(bool(settings),"settings write failed");
 auto manifest=output(stage/"simulation.json");manifest<<"{\n  \"format\": \"gsim.simulation/v1\",\n  \"seed\": "<<seed<<",\n  \"animals\": "<<n<<",\n  \"genotyped\": "<<ng<<",\n  \"markers\": "<<m<<",\n  \"traits\": "<<q<<",\n  \"observation_units\": "<<blocks<<",\n  \"observed_scalars\": "<<scalars<<",\n  \"records_per_individual\": "<<visits<<",\n  \"visits\": "<<visits<<",\n  \"missing_rate\": "<<missing<<",\n  \"maximum_memory_bytes\": "<<budget<<",\n  \"rng\": \"existing identity-keyed SplitMix64 streams and meiosis event sampler\",\n  \"qualification\": \"Frozen pedigree/founder/effect workload; correlated Gaussian residuals. Marker-derived truth, not an exact Gaussian H draw or variance-recovery qualification.\",\n  \"files_fnv1a64\": {";
 bool first=true;for(auto name:{"genotypes.bed","genotypes.bim","genotypes.fam","pedigree.csv","true_breeding_values.csv","phenotypes.csv","record_truth.csv","effective_settings.json"}){if(!first)manifest<<',';first=false;manifest<<'\n'<<"    "<<std::quoted(name)<<": "<<std::quoted(std::to_string(checksum(stage/name)));}manifest<<"\n  },\n  \"inputs_fnv1a64\": {";first=true;for(auto k:{"pedigree","variants","founders","effects"}){if(!first)manifest<<',';first=false;manifest<<'\n'<<"    "<<std::quoted(k)<<": "<<std::quoted(std::to_string(checksum(c.at(k))));}manifest<<"\n  }\n}\n";manifest.close();need(bool(manifest),"manifest write failed");for(auto k:{"pedigree","variants","founders","effects"})need(input_identity.at(k)==checksum(c.at(k)),"simulation input changed during reading");fs::rename(stage,out);
}
}
