#include "bed_reader.h"
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace breeding {
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
std::string string(SEXP x,R_xlen_t i) {
  require(TYPEOF(x)==STRSXP && i<XLENGTH(x) && STRING_ELT(x,i)!=NA_STRING,"invalid string identity");
  return Rf_translateCharUTF8(STRING_ELT(x,i));
}
std::uint64_t mix(std::uint64_t x) {
  x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);return x^(x>>31);
}
std::uint64_t hash(const std::string& s) {
  std::uint64_t h=UINT64_C(14695981039346656037);for(unsigned char c:s){h^=c;h*=UINT64_C(1099511628211);}return h;
}
std::uint64_t seed(SEXP x) {
  require(TYPEOF(x)==REALSXP && XLENGTH(x)==1 && R_FINITE(REAL(x)[0]) && REAL(x)[0]>=0 &&
    REAL(x)[0]<=9007199254740991. && REAL(x)[0]==std::floor(REAL(x)[0]),"invalid identity-stream seed");
  return static_cast<std::uint64_t>(REAL(x)[0]);
}
struct Rng {
  using result_type=std::uint64_t;
  std::uint64_t state;
  static constexpr result_type min(){return 0;} static constexpr result_type max(){return UINT64_MAX;}
  result_type operator()(){state+=UINT64_C(0x9e3779b97f4a7c15);return mix(state);}
  double uniform(){return (static_cast<double>((*this)()>>11)+.5)/9007199254740992.;}
  double normal(){return std::sqrt(-2*std::log(uniform()))*std::cos(6.2831853071795864769*uniform());}
};
Rng stream(std::uint64_t s,const std::string& id,const std::string& component,const std::string& domain) {
  return {mix(mix(s)^mix(hash(id))^mix(hash(component)^UINT64_C(0x13198a2e03707344))^
    mix(hash(domain)^UINT64_C(0xa4093822299f31d0)))};
}
int integer(SEXP x) {require(TYPEOF(x)==INTSXP && XLENGTH(x)==1 && INTEGER(x)[0]!=NA_INTEGER,"invalid integer scalar");return INTEGER(x)[0];}
void matrix(SEXP x,int rows,int columns) {
  require(TYPEOF(x)==REALSXP && Rf_isMatrix(x),"expected numeric matrix");
  SEXP d=Rf_getAttrib(x,R_DimSymbol);require(INTEGER(d)[0]==rows && INTEGER(d)[1]==columns,"matrix dimensions do not align");
  for(R_xlen_t i=0;i<XLENGTH(x);++i)require(R_FINITE(REAL(x)[i]),"nonfinite matrix value");
}
SEXP named_list(const std::vector<const char*>& names) {
  SEXP out=PROTECT(Rf_allocVector(VECSXP,names.size())),labels=PROTECT(Rf_allocVector(STRSXP,names.size()));
  for(std::size_t i=0;i<names.size();++i)SET_STRING_ELT(labels,i,Rf_mkChar(names[i]));
  Rf_setAttrib(out,R_NamesSymbol,labels);UNPROTECT(2);return out;
}
}

extern "C" SEXP C_gsim_breeding_checksum(SEXP bytes) {
  if(TYPEOF(bytes)!=RAWSXP)Rf_error("native breeding checksum requires serialized bytes");
  std::uint64_t value=104729;
  for(R_xlen_t i=0;i<XLENGTH(bytes);++i)value=(value*131+RAW(bytes)[i]+1)%UINT64_C(2147483629);
  char result[9];std::snprintf(result,sizeof(result),"%08x",static_cast<unsigned>(value));
  return Rf_mkString(result);
}

extern "C" SEXP C_gsim_breeding_order(SEXP animal,SEXP sire,SEXP dam) {
  try {
    using namespace breeding;const int n=static_cast<int>(XLENGTH(animal));
    require(TYPEOF(animal)==STRSXP && n>0 && TYPEOF(sire)==STRSXP && TYPEOF(dam)==STRSXP &&
      XLENGTH(sire)==n && XLENGTH(dam)==n,"pedigree identifiers do not align");
    std::unordered_map<std::string,int> index;for(int i=0;i<n;++i) {
      auto id=string(animal,i);require(!id.empty() && index.emplace(id,i).second,"animal IDs must be unique and nonempty");
    }
    std::vector<std::vector<int>> children(n);std::vector<int> degree(n),generation(n,1),order;
    for(int i=0;i<n;++i)for(SEXP side:{sire,dam})if(STRING_ELT(side,i)!=NA_STRING) {
      auto at=index.find(string(side,i));require(at!=index.end(),"parent ID is absent from supplied table");
      require(at->second!=i,"an animal cannot be its own parent");children[at->second].push_back(i);++degree[i];
    }
    for(int i=0;i<n;++i)if(!degree[i])order.push_back(i);
    for(std::size_t at=0;at<order.size();++at)for(int child:children[order[at]]) {
      generation[child]=std::max(generation[child],generation[order[at]]+1);
      if(--degree[child]==0)order.push_back(child);
    }
    require(order.size()==static_cast<std::size_t>(n),"pedigree contains a parentage cycle");
    SEXP out=PROTECT(named_list({"order","generation"}));
    SEXP a=PROTECT(Rf_allocVector(INTSXP,n)),b=PROTECT(Rf_allocVector(INTSXP,n));
    for(int i=0;i<n;++i){INTEGER(a)[i]=order[i]+1;INTEGER(b)[i]=generation[i];}
    SET_VECTOR_ELT(out,0,a);SET_VECTOR_ELT(out,1,b);UNPROTECT(3);return out;
  }catch(const std::exception& e){Rf_error("native pedigree construction: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_select(SEXP ids,SEXP strata,SEXP quotas,SEXP scores,SEXP random,SEXP sd) {
  try {
    using namespace breeding;const auto s=seed(sd);const int n=XLENGTH(ids),groups=XLENGTH(quotas);
    require(n>0 && TYPEOF(strata)==INTSXP && XLENGTH(strata)==n && TYPEOF(quotas)==INTSXP &&
      TYPEOF(scores)==REALSXP && XLENGTH(scores)==n,"selection vectors do not align");
    const bool shuffle=integer(random)!=0;std::vector<std::vector<int>> by(groups);
    std::vector<std::uint64_t> keys(n);std::vector<std::string> names(n);std::unordered_set<std::string> unique;
    for(int i=0;i<n;++i) {
      names[i]=string(ids,i);require(unique.insert(names[i]).second && !names[i].empty(),"duplicate or empty selection IDs");
      const int h=INTEGER(strata)[i]-1;require(h>=0 && h<groups,"invalid selection stratum");
      require(shuffle || R_FINITE(REAL(scores)[i]),"nonfinite selection score");by[h].push_back(i);
      keys[i]=stream(s,names[i],"rank","sampling")();
    }
    std::vector<int> selected;
    for(int h=0;h<groups;++h) {
      const int k=INTEGER(quotas)[h];require(k!=NA_INTEGER && k>=0 && k<=static_cast<int>(by[h].size()),"quota exceeds eligible stratum");
      auto compare=[&](int a,int b) {if(shuffle && keys[a]!=keys[b])return keys[a]<keys[b];
        if(!shuffle && REAL(scores)[a]!=REAL(scores)[b])return REAL(scores)[a]>REAL(scores)[b];
        return names[a]<names[b];};
      std::partial_sort(by[h].begin(),by[h].begin()+k,by[h].end(),compare);
      selected.insert(selected.end(),by[h].begin(),by[h].begin()+k);
    }
    SEXP out=PROTECT(Rf_allocVector(INTSXP,selected.size()));
    for(std::size_t i=0;i<selected.size();++i)INTEGER(out)[i]=selected[i]+1;
    UNPROTECT(1);return out;
  }catch(const std::exception& e){Rf_error("native selection: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_mate(SEXP child,SEXP sires,SEXP dams,SEXP sd) {
  try {
    using namespace breeding;const auto s=seed(sd);const int n=XLENGTH(child);
    require(n>0 && XLENGTH(sires)>0 && XLENGTH(dams)>0,"mating requires offspring and both parent pools");
    std::vector<std::string> father,mother;for(R_xlen_t i=0;i<XLENGTH(sires);++i)father.push_back(string(sires,i));
    for(R_xlen_t i=0;i<XLENGTH(dams);++i)mother.push_back(string(dams,i));
    std::sort(father.begin(),father.end());std::sort(mother.begin(),mother.end());
    SEXP out=PROTECT(named_list({"sire","dam","sex"}));
    for(int j=0;j<3;++j){SEXP x=PROTECT(Rf_allocVector(STRSXP,n));SET_VECTOR_ELT(out,j,x);UNPROTECT(1);}
    for(int i=0;i<n;++i){auto id=string(child,i);auto a=stream(s,id,"sire","mating"),b=stream(s,id,"dam","mating"),c=stream(s,id,"sex","mating");
      // Rejection sampling avoids modulo bias for unequal parent-pool sizes.
      auto draw=[](Rng& r,std::uint64_t count){const auto limit=-count%count;std::uint64_t x;do{x=r();}while(x<limit);return x%count;};
      SET_STRING_ELT(VECTOR_ELT(out,0),i,Rf_mkCharCE(father[draw(a,father.size())].c_str(),CE_UTF8));
      SET_STRING_ELT(VECTOR_ELT(out,1),i,Rf_mkCharCE(mother[draw(b,mother.size())].c_str(),CE_UTF8));
      SET_STRING_ELT(VECTOR_ELT(out,2),i,Rf_mkChar(c.uniform()<.5?"M":"F"));
    }UNPROTECT(1);return out;
  }catch(const std::exception& e){Rf_error("native mating: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_records(SEXP genetic,SEXP owner,SEXP design,SEXP fixed,SEXP animal,
  SEXP unit,SEXP trait,SEXP trait_names,SEXP residual_factor,SEXP permanent_factor,SEXP sd) {
  try {
    using namespace breeding;const auto s=seed(sd);const int n=XLENGTH(owner),t=XLENGTH(trait_names);
    require(n>0 && t>0 && t<=64 && TYPEOF(owner)==INTSXP && TYPEOF(trait)==INTSXP && XLENGTH(trait)==n &&
      TYPEOF(fixed)==REALSXP && XLENGTH(fixed)==n && XLENGTH(animal)==n && XLENGTH(unit)==n,"record inputs do not align");
    require(Rf_isMatrix(genetic) && TYPEOF(genetic)==REALSXP,"genomic values must be a numeric matrix");
    SEXP dims=Rf_getAttrib(genetic,R_DimSymbol);const int animals=INTEGER(dims)[0],q=INTEGER(dims)[1];
    matrix(genetic,animals,q);matrix(design,n,q);matrix(residual_factor,t,t);matrix(permanent_factor,t,t);
    SEXP out=PROTECT(Rf_allocMatrix(REALSXP,n,4));
    for(int i=0;i<n;++i) {
      require(INTEGER(owner)[i]>=1 && INTEGER(owner)[i]<=animals && INTEGER(trait)[i]>=1 && INTEGER(trait)[i]<=t,"record owner/trait out of range");
      const int a=INTEGER(owner)[i]-1,h=INTEGER(trait)[i]-1;
      require(R_FINITE(REAL(fixed)[i]),"nonfinite fixed contribution");double g=0,r=0,p=0;
      for(int j=0;j<q;++j)g+=REAL(genetic)[a+j*animals]*REAL(design)[i+j*n];
      for(int j=0;j<t;++j) {
        auto e=stream(s,string(unit,i),string(trait_names,j),"residual"),v=stream(s,string(animal,i),string(trait_names,j),"permanent");
        r+=REAL(residual_factor)[h+j*t]*e.normal();p+=REAL(permanent_factor)[h+j*t]*v.normal();
      }
      REAL(out)[i]=g;REAL(out)[i+n]=p;REAL(out)[i+2*n]=r;REAL(out)[i+3*n]=g+p+r+REAL(fixed)[i];
    }UNPROTECT(1);return out;
  }catch(const std::exception& e){Rf_error("native genomic records: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_normalize(SEXP group,SEXP weights,SEXP groups) {
  try {
    using namespace breeding;const int p=integer(groups),n=XLENGTH(group);require(p>0 && TYPEOF(group)==INTSXP &&
      TYPEOF(weights)==REALSXP && XLENGTH(weights)==n,"pool weights do not align");std::vector<double> sum(p);
    for(int i=0;i<n;++i){const int h=INTEGER(group)[i]-1;const double w=REAL(weights)[i];require(h>=0 && h<p && R_FINITE(w) && w>=0,"invalid pool weight");sum[h]+=w;}
    for(double x:sum)require(R_FINITE(x) && x>0,"every pool must have positive total contribution");
    SEXP out=PROTECT(Rf_allocVector(REALSXP,n));for(int i=0;i<n;++i)REAL(out)[i]=REAL(weights)[i]/sum[INTEGER(group)[i]-1];UNPROTECT(1);return out;
  }catch(const std::exception& e){Rf_error("native pool normalization: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_pool_bed(SEXP file,SEXP samples,SEXP markers,SEXP rows,SEXP groups,
  SEXP weights,SEXP pool_ids,SEXP columns,SEXP marker_ids,SEXP depth,SEXP sd) {
  try {
    using namespace breeding;const auto s=seed(sd);const int n=integer(samples),m=integer(markers),p=XLENGTH(pool_ids),c=XLENGTH(columns),k=XLENGTH(rows),d=integer(depth);
    require(n>0 && m>0 && p>0 && c>0 && c<=64 && d>=0 && TYPEOF(rows)==INTSXP && TYPEOF(groups)==INTSXP &&
      TYPEOF(columns)==INTSXP && TYPEOF(weights)==REALSXP && XLENGTH(groups)==k && XLENGTH(weights)==k && XLENGTH(marker_ids)==c,"invalid pooled BED dimensions");
    for(int i=0;i<k;++i)require(INTEGER(rows)[i]>=1 && INTEGER(rows)[i]<=n && INTEGER(groups)[i]>=1 && INTEGER(groups)[i]<=p &&
      R_FINITE(REAL(weights)[i]) && REAL(weights)[i]>=0,"invalid pooled BED contribution");
    gsim::native::BedReader reader(string(file,0),n,m);SEXP out=PROTECT(Rf_allocMatrix(REALSXP,p,c));std::fill(REAL(out),REAL(out)+XLENGTH(out),0.);
    SEXP variance=PROTECT(Rf_allocMatrix(REALSXP,p,c));std::fill(REAL(variance),REAL(variance)+XLENGTH(variance),0.);
    for(int j=0;j<c;++j){const int col=INTEGER(columns)[j]-1;require(col>=0 && col<m,"pooled marker out of range");reader.read_record(col);
      for(int i=0;i<k;++i)if(REAL(weights)[i]>0){const int dosage=reader.dosage(INTEGER(rows)[i]-1);require(dosage>=0,"missing genotype in pooled DNA contribution");REAL(out)[INTEGER(groups)[i]-1+j*p]+=.5*REAL(weights)[i]*dosage;}
      for(int h=0;h<p;++h){double& frequency=REAL(out)[h+j*p];require(frequency>=-1e-12 && frequency<=1+1e-12,"pool weights are not normalized");frequency=std::min(1.,std::max(0.,frequency));
        if(d){REAL(variance)[h+j*p]=frequency*(1-frequency)/d;
          auto rng=stream(s,string(pool_ids,h),string(marker_ids,j),"assay");std::binomial_distribution<int> binomial(d,frequency);frequency=static_cast<double>(binomial(rng))/d;}}
    }Rf_setAttrib(out,Rf_install("assay_variance"),variance);UNPROTECT(2);return out;
  }catch(const std::exception& e){Rf_error("native pooled BED: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_pool_covariance(SEXP pool,SEXP weight,SEXP animal,SEXP unit,SEXP trait,
  SEXP residual_factor,SEXP permanent_factor,SEXP groups,SEXP left,SEXP right,SEXP diagonal) {
  try {
    using namespace breeding;const int p=integer(groups),n=XLENGTH(pool);require(p>0 && TYPEOF(pool)==INTSXP && TYPEOF(trait)==INTSXP &&
      TYPEOF(weight)==REALSXP && XLENGTH(weight)==n && XLENGTH(animal)==n && XLENGTH(unit)==n && XLENGTH(trait)==n,"pool noise inputs do not align");
    require(Rf_isMatrix(residual_factor),"expected noise factor");const int t=INTEGER(Rf_getAttrib(residual_factor,R_DimSymbol))[0];
    require(t>0 && t<=64,"noise traits exceed supported bound");matrix(residual_factor,t,t);matrix(permanent_factor,t,t);
    using Terms=std::unordered_map<std::string,std::vector<double>>;
    std::vector<Terms> residual(p),permanent(p);
    for(int i=0;i<n;++i) {
      const int h=INTEGER(pool)[i]-1,a=INTEGER(trait)[i]-1;const double w=REAL(weight)[i];require(h>=0 && h<p && a>=0 && a<t && R_FINITE(w),"invalid noise contribution");
      auto accumulate=[&](Terms& terms,const std::string& id,SEXP factor){auto& v=terms[id];if(v.empty())v.resize(t);for(int j=0;j<t;++j)v[j]+=w*REAL(factor)[a+j*t];};
      accumulate(residual[h],string(unit,i),residual_factor);accumulate(permanent[h],string(animal,i),permanent_factor);
    }
    auto product=[&](const Terms& a,const Terms& b){double value=0;for(const auto& term:a){auto found=b.find(term.first);if(found!=b.end())for(int j=0;j<t;++j)value+=term.second[j]*found->second[j];}return value;};
    const bool diag=integer(diagonal)!=0;const int l=XLENGTH(left),r=XLENGTH(right);
    require(TYPEOF(left)==INTSXP && TYPEOF(right)==INTSXP && (diag || (l<=256 && r<=256)),"covariance query exceeds 256 by 256");
    SEXP out=PROTECT(diag?Rf_allocVector(REALSXP,l):Rf_allocMatrix(REALSXP,l,r));
    for(int i=0;i<l;++i){int a=INTEGER(left)[i]-1;require(a>=0 && a<p,"covariance pool index out of range");
      for(int j=0;j<(diag?1:r);++j){int b=diag?a:INTEGER(right)[j]-1;require(b>=0 && b<p,"covariance pool index out of range");
        REAL(out)[i+(diag?0:j*l)]=product(residual[a],residual[b])+product(permanent[a],permanent[b]);}}
    UNPROTECT(1);return out;
  }catch(const std::exception& e){Rf_error("native pooled covariance: %s",e.what());}return R_NilValue;
}

extern "C" SEXP C_gsim_breeding_pool_values(SEXP pool,SEXP weight,SEXP value,SEXP groups) {
  try {
    using namespace breeding;const int p=integer(groups),n=XLENGTH(pool);
    require(p>0 && TYPEOF(pool)==INTSXP && TYPEOF(weight)==REALSXP && TYPEOF(value)==REALSXP &&
      XLENGTH(weight)==n && XLENGTH(value)==n,"pool values do not align");
    SEXP out=PROTECT(Rf_allocVector(REALSXP,p));std::fill(REAL(out),REAL(out)+p,0.);
    for(int i=0;i<n;++i){const int h=INTEGER(pool)[i];const double w=REAL(weight)[i],v=REAL(value)[i];
      require(h>=1 && h<=p && R_FINITE(w) && w>=0 && R_FINITE(v),"invalid pooled value contribution");REAL(out)[h-1]+=w*v;}
    UNPROTECT(1);return out;
  }catch(const std::exception& e){Rf_error("native pooled values: %s",e.what());}return R_NilValue;
}
