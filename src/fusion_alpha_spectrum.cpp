#include "fusion_alpha_spectrum.h"
#include "fusion_alpha_amplitudes.h"
#include "fusion_products.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
namespace {
constexpr double pi=3.1415926535897932384626433832795,mev=1.602176634e-13;
using C=std::complex<double>;
struct Node {double x,w;};
std::vector<Node> gauss(int n){
 std::vector<Node> a(n);
 for(int i=0;i<(n+1)/2;++i){double x=std::cos(pi*(i+.75)/(n+.5)),derivative=0;
  bool converged=false;
  for(int it=0;it<64;++it){double p=1,prev=0;
   for(int j=1;j<=n;++j){double next=((2*j-1)*x*p-(j-1)*prev)/j;prev=p;p=next;}
   derivative=n*(x*p-prev)/(x*x-1);double change=p/derivative;x-=change;
   if(std::abs(change)<4e-16){converged=true;break;}
  }
  if(!converged || !std::isfinite(x) || std::abs(x)>=1)
   throw std::runtime_error("Gauss-Legendre root did not converge");
  // Recompute derivative at the final abscissa.
  double p=1,prev=0;for(int j=1;j<=n;++j){double next=((2*j-1)*x*p-(j-1)*prev)/j;prev=p;p=next;}
  derivative=n*(x*p-prev)/(x*x-1);double w=2/((1-x*x)*derivative*derivative);
  if(!std::isfinite(w) || w<=0)throw std::runtime_error("Invalid quadrature weight");
  a[i]={-x,w};a[n-1-i]={x,w};
 }
 long double weight=0;for(const auto &v:a)weight+=v.w;
 if(std::abs(weight-2)>1e-12)throw std::runtime_error("Quadrature weight closure");
 return a;
}
struct Event {std::array<C,5> a{},b{};double e[3],measure;};
double norm(const std::array<C,5>&a){double s=0;for(C x:a)s+=std::norm(x);return s;}
}
static int spectrum_impl(int mode,int policy,double A,double cutoff,
 double k,double phase,int nq,int nc,int cells,const double *edges,double *birth,
 fusion_alpha_spectrum_v1 *out){
 if(out)*out={};
 if(cells>0 && birth)std::fill(birth,birth+cells,0.);
 if(!birth || !out)return PB11_STATUS_NULL_OUTPUT;
 if(!edges || cells<1 || nq<4 || nc<4 || nq>1024 || nc>1024 ||
  (policy!=0 && policy!=1) || (mode!=1 && mode!=2 && mode!=3 && mode!=13) || !std::isfinite(A) ||
  !std::isfinite(k) || !std::isfinite(phase) || !std::isfinite(cutoff))return PB11_STATUS_INVALID_ARGUMENT;
 if(A<=0 || A>12*mev || k<0 || k>1 || cutoff<.001*mev || cutoff>.01*mev)return PB11_STATUS_OUT_OF_RANGE;
 for(int i=0;i<=cells;++i)
  if(!std::isfinite(edges[i]) || edges[i]<0 || (i>0 && edges[i]<=edges[i-1]))
   return PB11_STATUS_INVALID_ARGUMENT;
 try{
  const auto qnodes=gauss(nq),cnodes=gauss(nc);std::vector<Event> events;events.reserve(nq*nc);
  long double n1=0,n3=0;int pruned_events=0;
  for(const auto &u:qnodes){double theta=(u.x+1)*pi/4,s=std::sin(theta),q=A*s*s;
   double dq=A*std::sin(2*theta)*u.w*pi/4;
   for(const auto &v:cnodes){Event e{};fusion_alpha_amplitudes_v1 a{},b{};int pa=0,pb=0;
    int status=fusion_c_alpha_amplitudes_fsci_cutoff(mode==13?1:mode,policy,A,q,v.x,cutoff,&a,&pa);if(status)return status;
    if(mode==13){status=fusion_c_alpha_amplitudes_fsci_cutoff(3,policy,A,q,v.x,cutoff,&b,&pb);if(status)return status;}
    e.measure=a.phase_space_J*dq*v.w;
    for(int M=0;M<5;++M){e.a[M]={a.sym_real[M],a.sym_imag[M]};e.b[M]={b.sym_real[M],b.sym_imag[M]};}
    n1+=static_cast<long double>(e.measure)*norm(e.a);n3+=static_cast<long double>(e.measure)*norm(e.b);
    if(pa || pb)++pruned_events;
    // Positive squared nonrelativistic momenta avoid endpoint cancellation.
    double p0=std::sqrt(4*(A-q)/3),star=std::sqrt(q),px=star*std::sqrt((1-v.x)*(1+v.x));
    e.e[0]=p0*p0/2;e.e[1]=(px*px+std::pow(-p0/2-star*v.x,2))/2;e.e[2]=(px*px+std::pow(-p0/2+star*v.x,2))/2;
    events.push_back(e);
   }
  }
  if(!(n1>0) || (mode==13 && !(n3>0)))return PB11_STATUS_NUMERICAL_FAILURE;
  const C factor=mode==13?std::sqrt((1-k)*static_cast<double>(n1/n3))*std::polar(1.,phase):C(0,0);
  std::vector<double> rates(events.size());long double total=0;
  for(std::size_t i=0;i<events.size();++i){std::array<C,5> combined=events[i].a;
   if(mode==13)for(int M=0;M<5;++M)combined[M]=std::sqrt(k)*events[i].a[M]+factor*events[i].b[M];
   rates[i]=events[i].measure*norm(combined);total+=rates[i];}
  if(!(total>0) || !std::isfinite(total))return PB11_STATUS_NUMERICAL_FAILURE;
  std::vector<double> packet_e,packet_w;packet_e.reserve(3*events.size());packet_w.reserve(3*events.size());
  for(std::size_t i=0;i<events.size();++i)for(double e:events[i].e){packet_e.push_back(e);packet_w.push_back(static_cast<double>(rates[i]/total));}
  std::vector<double> mapped(cells);fusion_birth_mapping_v1 ledger{};
  // Unit reaction rate1 m^-3 s^-1; divide all outputs by that unit rate.
  int status=fusion_c_map_birth_packets(cells,edges,static_cast<int>(packet_e.size()),packet_e.data(),packet_w.data(),mapped.data(),&ledger);if(status)return status;
  fusion_alpha_spectrum_v1 result{};
  result.mapped_number=ledger.mapped_number_m3_s;result.mapped_energy_J=ledger.mapped_energy_W_m3;
  result.below_number=ledger.below_number_m3_s;result.below_energy_J=ledger.below_energy_W_m3;
  result.above_number=ledger.above_number_m3_s;result.above_energy_J=ledger.above_energy_W_m3;
  result.number_residual=ledger.input_number_m3_s-3;result.energy_residual_J=ledger.input_energy_W_m3-A;
  result.normalization_J2=static_cast<double>(total);
  result.l1_normalization_J2=mode==13?static_cast<double>(n1):0;result.l3_normalization_J2=mode==13?static_cast<double>(n3):0;
  result.quadrature_events=nq*nc;result.pruned_events=pruned_events;
  if(std::abs(result.number_residual)>3e-12 || std::abs(result.energy_residual_J)>1e-12*A ||
   !std::isfinite(result.normalization_J2) || result.normalization_J2==0)return PB11_STATUS_NUMERICAL_FAILURE;
  std::copy(mapped.begin(),mapped.end(),birth);*out=result;return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}

extern "C" int fusion_c_alpha_spectrum_grid(int mode,double A,double cutoff,
 double k,double phase,int nq,int nc,int cells,const double *edges,double *birth,
 fusion_alpha_spectrum_v1 *out){
 return spectrum_impl(mode,0,A,cutoff,k,phase,nq,nc,cells,edges,birth,out);
}
extern "C" int fusion_c_alpha_spectrum_model_grid(int mode,int policy,double A,
 double cutoff,double k,double phase,int nq,int nc,int cells,const double *edges,
 double *birth,fusion_alpha_spectrum_v1 *out){
 return spectrum_impl(mode,policy,A,cutoff,k,phase,nq,nc,cells,edges,birth,out);
}
