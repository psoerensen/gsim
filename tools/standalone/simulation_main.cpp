#include "simulation_workflow.h"
#include <fstream>
#include <iostream>
#include <iterator>
int main(int argc,char**argv){try{
 if(argc==2&&std::string(argv[1])=="--version"){std::cout<<"gsim standalone simulation contract 1\n";return 0;}
 if(argc!=4||std::string(argv[1])!="--wire")throw std::runtime_error("usage: gsim_native --wire recipe.txt output-directory");
 std::ifstream in(argv[2],std::ios::binary);if(!in)throw std::runtime_error("cannot open simulation recipe");
 gsim::native::simulate_recipe(std::string(std::istreambuf_iterator<char>(in),{}),argv[3]);return 0;
}catch(const std::exception&e){std::cerr<<"gsim: "<<e.what()<<'\n';return 1;}}
