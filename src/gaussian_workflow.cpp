#include "gaussian_draws.h"
#include "breeding_rng.h"
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <limits>
#include <chrono>
#include <thread>
namespace gsim::native {namespace {
using U=std::uint64_t;using V=std::vector<double>;namespace fs=std::filesystem;
void need(bool ok,const std::string& text){if(!ok)throw std::runtime_error(text);}
std::map<std::string,std::string> parse(std::istream& in,const std::string& magic){std::string line;std::getline(in,line);if(!line.empty()&&line.back()=='\r')line.pop_back();need(line==magic,"unsupported Gaussian simulation specification");std::map<std::string,std::string> out;std::string k,v;while(in>>k){need(bool(in>>std::quoted(v)),"truncated Gaussian simulation setting");need(out.emplace(k,v).second,"duplicate Gaussian simulation setting");}need(in.eof(),"invalid Gaussian simulation setting");return out;}
U integer(const std::string& s){need(!s.empty()&&s.find_first_not_of("0123456789")==s.npos,"Gaussian integer must contain decimal digits");std::size_t at;auto value=std::stoull(s,&at);need(at==s.size(),"invalid Gaussian integer");return value;}
U hash(const fs::path& p){std::ifstream f(p,std::ios::binary);need(bool(f),"cannot checksum Gaussian input");U h=14695981039346656037ULL;char b[65536];while(f){f.read(b,sizeof b);for(std::streamsize i=0;i<f.gcount();++i){h^=static_cast<unsigned char>(b[i]);h*=1099511628211ULL;}}need(f.eof(),"Gaussian checksum read error");return h;}
U get(std::istream& f){U value=0;for(U i=0;i<8;++i){int c=f.get();need(c!=EOF,"truncated Gaussian binary array");value|=U(static_cast<unsigned char>(c))<<(8*i);}return value;}
V array(const fs::path& p,U count){need(count<=(std::numeric_limits<U>::max()-16)/8&&fs::file_size(p)==16+8*count,"Gaussian binary array file length differs");std::ifstream f(p,std::ios::binary);char magic[8];need(bool(f.read(magic,8))&&std::memcmp(magic,"GSVEC001",8)==0&&get(f)==count,"invalid Gaussian binary header");V out(count);for(auto& v:out){auto bits=get(f);std::memcpy(&v,&bits,8);need(std::isfinite(v),"nonfinite Gaussian input array");}return out;}
std::ofstream output(const fs::path& p){std::ofstream f(p,std::ios::binary);f.imbue(std::locale::classic());f<<std::setprecision(17);need(bool(f),"cannot open Gaussian output");return f;}
std::string csv(const std::string& s){std::string out="\"";for(char c:s){out+=c;if(c=='"')out+='"';}return out+'"';}
void normals(U n,U seed,const char* domain,V& z){z.resize(n);for(U i=0;i<n;++i)z[i]=breeding::stream(seed,domain,std::to_string(i),"gaussian_model").normal();}
// Windows scanners can briefly retain a newly closed file handle. Keep
// publication atomic and bounded; never overwrite or copy a partial dataset.
void publish(const fs::path& stage,const fs::path& out){
 for(unsigned attempt=0;;++attempt){
  need(!fs::exists(out),"preserve existing Gaussian simulation output");
  std::error_code error;fs::rename(stage,out,error);if(!error)return;
#ifdef _WIN32
  if(attempt<19 && (error==std::errc::permission_denied || error==std::errc::device_or_resource_busy)){
   std::this_thread::sleep_for(std::chrono::milliseconds(25));continue;
  }
#endif
  throw fs::filesystem_error("cannot publish Gaussian simulation",stage,out,error);
 }
}
V apply(const V& a,U n,U m,const V& x){need(a.size()==n*m&&x.size()==m,"Gaussian product dimension mismatch");V y(n,0.);for(U i=0;i<n;++i)for(U j=0;j<m;++j)y[i]+=a[i+j*n]*x[j];for(double v:y)need(std::isfinite(v),"nonfinite Gaussian product");return y;}
}
void simulate_gaussian_recipe(const std::string& wire,const std::string& destination){
 need(wire.size()<=1024*1024,"Gaussian control input exceeds metadata allowance");
 std::istringstream top(wire);auto controls=parse(top,"GSIM_GAUSSIAN 1");for(const auto& kv:controls)need(kv.first=="specification"||kv.first=="seed"||kv.first=="maximum_memory_bytes","unknown Gaussian simulation control");
 need(controls.count("specification"),"missing prepared Gaussian specification");const U seed=integer(controls.count("seed")?controls.at("seed"):"1"),budget=integer(controls.count("maximum_memory_bytes")?controls.at("maximum_memory_bytes"):"4294967296");need(seed<=9007199254740989ULL,"Gaussian simulation seed exceeds exact range");
 fs::path specification=fs::absolute(controls.at("specification")).lexically_normal();need(fs::file_size(specification)<=1024*1024,"Gaussian specification exceeds metadata allowance");std::ifstream input(specification);need(bool(input),"cannot open Gaussian specification");auto settings=parse(input,"GSIM_GAUSSIAN_INPUTS 1");
 const std::set<std::string> keys{"latent_dimension","coefficient_dimension","primary_dimension","record_count","prior_covariance","coefficient_projection","coefficient_mean","record_design","residual_covariance","record_labels","provenance"};for(const auto& kv:settings)need(keys.count(kv.first),"unknown Gaussian specification key");for(const auto& key:keys)need(settings.count(key),"incomplete Gaussian specification");
 const U latent=integer(settings.at("latent_dimension")),d=integer(settings.at("coefficient_dimension")),primary=integer(settings.at("primary_dimension")),n=integer(settings.at("record_count"));
 need(latent<=1024&&d>0&&d<=4096&&primary>0&&primary<=d&&n>0&&n<=1024,"Gaussian simulation exceeds bounded model dimensions");
 const long double estimate=32.L*latent*latent+32.L*n*n+16.L*d*latent+16.L*n*d+64.L*(n+d+latent);need(estimate<=budget,"Gaussian simulation working-memory admission exceeds maximum_memory_bytes");
 std::map<std::string,fs::path> files;std::map<std::string,U> before;
 for(const auto& key:{"prior_covariance","coefficient_projection","coefficient_mean","record_design","residual_covariance","record_labels","provenance"}){auto path=fs::u8path(settings.at(key));if(path.is_relative())path=specification.parent_path()/path;files[key]=fs::absolute(path).lexically_normal();before[key]=hash(files[key]);}
 need(fs::file_size(files.at("record_labels"))<=32ULL*1024*1024&&fs::file_size(files.at("provenance"))<=32ULL*1024*1024,"Gaussian labels/provenance exceed metadata allowance");
 auto prior=gaussian_lower(array(files.at("prior_covariance"),latent*latent),latent),residual=gaussian_lower(array(files.at("residual_covariance"),n*n),n);
 auto transform=array(files.at("coefficient_projection"),d*latent),mean=array(files.at("coefficient_mean"),d),design=array(files.at("record_design"),n*d);
 struct Label{std::string record,trait,individual;};std::vector<Label> labels(n);std::set<std::pair<std::string,std::string>> identities;std::set<std::string> individuals,traits,units;
 std::ifstream names(files.at("record_labels"));std::string header;std::getline(names,header);if(!header.empty()&&header.back()=='\r')header.pop_back();need(header=="GSIM_GAUSSIAN_LABELS 1","invalid Gaussian record-label header");for(auto& label:labels){need(bool(names>>std::quoted(label.record)>>std::quoted(label.trait)>>std::quoted(label.individual))&&!label.record.empty()&&!label.trait.empty()&&!label.individual.empty(),"invalid/truncated Gaussian record labels");need(identities.emplace(label.record,label.trait).second,"duplicate Gaussian record/trait identity");individuals.insert(label.individual);traits.insert(label.trait);units.insert(label.record);}names>>std::ws;need(names.eof(),"extra Gaussian record labels");
 V z,e;normals(latent,seed,"coefficient",z);normals(n,seed+1,"residual",e);auto random=apply(prior,latent,latent,z),coefficient=apply(transform,d,latent,random);for(U i=0;i<d;++i)coefficient[i]+=mean[i];auto noise=apply(residual,n,n,e),signal=apply(design,n,d,coefficient),fixed=apply(design,n,d,mean);
 fs::path out=fs::absolute(destination).lexically_normal(),stage=out.string()+".staging";need(!fs::exists(out)&&!fs::exists(stage),"preserve existing Gaussian simulation output/staging directory");need(fs::is_directory(out.parent_path()),"Gaussian output parent must exist");fs::create_directory(stage);
 auto values=output(stage/"coefficients.csv");values<<"index,value\n";for(U i=0;i<d;++i)values<<i<<','<<coefficient[i]<<'\n';values.close();need(bool(values),"Gaussian coefficient output failed");
 auto records=output(stage/"record_truth.csv"),phenotypes=output(stage/"phenotypes.csv");records<<"index,record_id,individual,trait,random_effects,fixed,residual,phenotype\n";phenotypes<<"record_id,individual,trait,value\n";
 for(U i=0;i<n;++i){const auto& label=labels[i];const double y=signal[i]+noise[i];need(std::isfinite(y),"nonfinite Gaussian phenotype");records<<i<<','<<csv(label.record)<<','<<csv(label.individual)<<','<<csv(label.trait)<<','<<signal[i]-fixed[i]<<','<<fixed[i]<<','<<noise[i]<<','<<y<<'\n';phenotypes<<csv(label.record)<<','<<csv(label.individual)<<','<<csv(label.trait)<<','<<y<<'\n';}records.close();phenotypes.close();need(bool(records)&&bool(phenotypes),"Gaussian record output failed");
 for(auto item:{std::make_pair("latent_innovations.csv",&z),std::make_pair("residual_innovations.csv",&e)}){auto f=output(stage/item.first);f<<"index,standard_normal\n";for(U i=0;i<item.second->size();++i)f<<i<<','<<(*item.second)[i]<<'\n';f.close();need(bool(f),"Gaussian innovation output failed");}
 fs::copy_file(files.at("provenance"),stage/"model_provenance.json");
 auto effective=output(stage/"effective_settings.json");effective<<"{\n  \"schema\": \"gsim.gaussian_model/v1\",\n  \"seed\": "<<seed<<",\n  \"latent_dimension\": "<<latent<<",\n  \"coefficient_dimension\": "<<d<<",\n  \"primary_dimension\": "<<primary<<",\n  \"maximum_memory_bytes\": "<<budget<<"\n}\n";effective.close();need(bool(effective),"Gaussian settings output failed");
 auto manifest=output(stage/"simulation.json");manifest<<"{\n  \"format\": \"gsim.gaussian_model/v1\",\n  \"animals\": "<<individuals.size()<<",\n  \"traits\": "<<traits.size()<<",\n  \"observation_units\": "<<units.size()<<",\n  \"observed_scalars\": "<<n<<",\n  \"seed\": "<<seed<<",\n  \"maximum_memory_bytes\": "<<budget<<",\n  \"qualification\": \"Gaussian draw conditional on the declared design, genotype/pedigree inputs and covariance state; no fitted response used; exact to the declared model including explicit approximations.\",\n  \"inputs_fnv1a64\": {";bool first=true;for(const auto& item:before){if(!first)manifest<<',';first=false;manifest<<'\n'<<"    "<<std::quoted(item.first)<<": "<<std::quoted(std::to_string(item.second));}manifest<<"\n  }\n}\n";manifest.close();need(bool(manifest),"Gaussian manifest output failed");
 for(const auto& item:files)need(hash(item.second)==before.at(item.first),"Gaussian input changed during simulation");publish(stage,out);
}
}
