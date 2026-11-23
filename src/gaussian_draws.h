#ifndef GSIM_GAUSSIAN_DRAWS_H
#define GSIM_GAUSSIAN_DRAWS_H
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
namespace gsim::native {
// The existing standalone workflow's exact Cholesky arithmetic, shared by
// marker residuals and covariance-matched Gaussian model generation.
inline std::vector<double> gaussian_lower(const std::vector<double>& a,std::uint64_t q){
 auto need=[](bool x,const char* message){if(!x)throw std::runtime_error(message);};
 need(a.size()==q*q,"residual covariance dimensions mismatch");std::vector<double> l(q*q);
 for(std::uint64_t i=0;i<q;++i)for(std::uint64_t j=0;j<=i;++j){need(a[i*q+j]==a[j*q+i],"residual covariance must be symmetric");double v=a[i*q+j];for(std::uint64_t k=0;k<j;++k)v-=l[i+k*q]*l[j+k*q];if(i==j){need(v>0,"residual covariance must be positive definite");l[i+j*q]=std::sqrt(v);}else l[i+j*q]=v/l[j+j*q];}return l;
}
void simulate_gaussian_recipe(const std::string&,const std::string&);
}
#endif
