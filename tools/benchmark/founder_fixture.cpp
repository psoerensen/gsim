// Benchmark-only founder approximation. Inheritance remains gsim's existing
// identity-keyed packed meiosis implementation. This is not a coalescent model.
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc,char**argv){try{
 if(argc!=5)throw std::runtime_error("founder_fixture map count seed output");
 std::ifstream map(argv[1]);std::string id;int bp,chip;double cm;std::vector<double> pos;
 while(map>>id>>bp>>cm>>chip)pos.push_back(cm);
 if(pos.empty())throw std::runtime_error("empty map");
 int n=std::stoi(argv[2]);std::mt19937_64 rng(std::stoull(argv[3]));
 const auto uniform=[&](){return (double(rng()>>11)+.5)/9007199254740992.;};
 std::vector<double> frequency(pos.size());std::vector<bool> flip(pos.size());
 int region=-1;double p=0;for(std::size_t j=0;j<pos.size();++j){int r=int(pos[j]);if(r!=region){region=r;p=.005+.495*std::pow(uniform(),1.8);}frequency[j]=p;flip[j]=uniform()<.5;}
 std::ofstream out(argv[4],std::ios::binary);for(int i=0;i<n;++i){out<<"a"<<(i+1);for(int h=0;h<2;++h){std::string bits(pos.size(),'0');bool state=false;for(std::size_t j=0;j<pos.size();++j){double refresh=j?-std::expm1(-(pos[j]-pos[j-1])/.2):1.;if(!j||int(pos[j])!=int(pos[j-1])||uniform()<refresh)state=uniform()<frequency[j];bits[j]='0'+(state!=flip[j]);}out<<' '<<bits;}out<<'\n';}
 if(!out)throw std::runtime_error("founder write failed");return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
