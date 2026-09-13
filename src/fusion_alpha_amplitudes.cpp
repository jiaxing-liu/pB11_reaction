#include "fusion_alpha_amplitudes.h"
#include "fusion_nuclear_coulomb.h"
#include <algorithm>
#include <cmath>
#include <complex>
namespace {
using C=std::complex<double>;
constexpr double pi=3.1415926535897932384626433832795,mev=1.602176634e-13;
double fact(int n){double r=1;for(int k=2;k<=n;++k)r*=k;return r;}
double cg(int j1,int m1,int j2,int m2,int J,int M){
 if(m1+m2!=M || std::abs(m1)>j1 || std::abs(m2)>j2 || std::abs(M)>J)return 0;
 double pre=std::sqrt((2*J+1)*fact(J+j1-j2)*fact(J-j1+j2)*fact(j1+j2-J)/fact(j1+j2+J+1)*fact(J+M)*fact(J-M)*fact(j1-m1)*fact(j1+m1)*fact(j2-m2)*fact(j2+m2));
 double s=0;
 for(int k=0;k<=j1+j2;++k){int d[]={k,j1+j2-J-k,j1-m1-k,j2+m2-k,J-j2+m1+k,J-j1-m2+k};
  if(*std::min_element(d,d+6)<0)continue;
  double denominator=1;for(int n:d)denominator*=fact(n);s+=(k%2?-1:1)/denominator;}
 return pre*s;
}
C spherical(int l,int m,const double *p){
 double norm=std::hypot(std::hypot(p[0],p[1]),p[2]);
 double x=p[2]/norm; x=std::max(-1.,std::min(1.,x));
 int a=std::abs(m);double P=1;
 for(int k=1;k<=a;++k)P*=-(2*k-1)*std::sqrt((1-x)*(1+x));
 if(l>a){double prev=P;P=x*(2*a+1)*P;
  for(int k=a+2;k<=l;++k){double next=((2*k-1)*x*P-(k+a-1)*prev)/(k-a);prev=P;P=next;}}
 C y=std::sqrt((2*l+1)/(4*pi)*fact(l-a)/fact(l+a))*P*std::polar(1.,a*std::atan2(p[1],p[0]));
 return m<0?(a%2?-1.:1.)*std::conj(y):y;
}
}
static int alpha_impl(int l,double A_in,double q_in,double cosine,
 fusion_alpha_amplitudes_v1 *out,double cutoff_MeV,int *pruned){
 if(pruned)*pruned=0;
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(l<1 || l>3 || !std::isfinite(A_in) || !std::isfinite(q_in) || !std::isfinite(cosine))return PB11_STATUS_INVALID_ARGUMENT;
 if(A_in<=0 || q_in<0 || q_in>A_in || (cutoff_MeV==0 && (q_in==0 || q_in==A_in)) || std::abs(cosine)>1)return PB11_STATUS_OUT_OF_RANGE;
 const double A=A_in/mev,q=q_in/mev,m=3727.3794118;
 // Momenta are MeV/c in the declared nonrelativistic convention.
 const double p0=std::sqrt(4*m*(A-q)/3),star=std::sqrt(m*q),px=star*std::sqrt((1-cosine)*(1+cosine));
 const double p[3][3]={{0,0,p0},{px,0,-p0/2-star*cosine},{-px,0,-p0/2+star*cosine}};
 fusion_nuclear_coulomb_v1 at_level{};
 int status=fusion_c_nuclear_coulomb(4,3.129*mev,&at_level);if(status)return status;
 C amplitudes[3][5]{}; int skipped=0;
 const C ipow[4]={{1,0},{0,1},{-1,0},{0,-1}};
 for(int i=0;i<3;++i){int j=(i+1)%3,k=(i+2)%3;double relative[3],e23=0,ep=0;
  for(int d=0;d<3;++d){relative[d]=(p[j][d]-p[k][d])/2;e23+=relative[d]*relative[d]/m;ep+=3*p[i][d]*p[i][d]/(4*m);}
  // Explicit caller-selected numerical approximation: retain the event and
  // all other permutations. Never discard a whole event for one small channel.
  if(cutoff_MeV>0 && (ep<cutoff_MeV || e23<cutoff_MeV)){++skipped;continue;}
  fusion_nuclear_coulomb_v1 primary{},secondary{};
  status=fusion_c_nuclear_coulomb(l-1,ep*mev,&primary);if(status)return status;
  status=fusion_c_nuclear_coulomb(4,e23*mev,&secondary);if(status)return status;
  const double g2=1.075,P=std::exp(secondary.log_penetrability);
  C denominator(3.129-e23-g2*(secondary.shift-at_level.shift),-g2*P);
  C radial=2*std::sqrt(g2/(primary.rho*secondary.rho))*std::exp(.5*(primary.log_penetrability+secondary.log_penetrability))*
   C(primary.phase_real,primary.phase_imag)*C(secondary.phase_real,secondary.phase_imag)*ipow[l]*ipow[2]/denominator;
  for(int M=-2;M<=2;++M)for(int b=-2;b<=2;++b){int a=M-b;if(std::abs(a)>l)continue;
   amplitudes[i][M+2]+=cg(2,b,l,a,2,M)*spherical(l,a,p[i])*spherical(2,b,relative)*radial;}
 }
 fusion_alpha_amplitudes_v1 result{};
 for(int M=0;M<5;++M){C s=amplitudes[0][M]+amplitudes[1][M]+amplitudes[2][M];
  if(!std::isfinite(s.real()) || !std::isfinite(s.imag()))return PB11_STATUS_NUMERICAL_FAILURE;
  result.unsym_real[M]=amplitudes[0][M].real();result.unsym_imag[M]=amplitudes[0][M].imag();
  result.sym_real[M]=s.real();result.sym_imag[M]=s.imag();}
 result.phase_space_J=std::sqrt(q_in)*std::sqrt(A_in-q_in);
 *out=result;if(pruned)*pruned=skipped;return PB11_STATUS_OK;
}

extern "C" int fusion_c_alpha_amplitudes(int l,double A,double q,double cosine,
 fusion_alpha_amplitudes_v1 *out){return alpha_impl(l,A,q,cosine,out,0,nullptr);}
extern "C" int fusion_c_alpha_amplitudes_cutoff(int l,double A,double q,double cosine,
 double cutoff_J,fusion_alpha_amplitudes_v1 *out,int *pruned){
 if(out)*out={};
 if(pruned)*pruned=0;
 if(!out || !pruned)return PB11_STATUS_NULL_OUTPUT;
 constexpr double minimum=.001*mev,maximum=.01*mev;
 if(!std::isfinite(cutoff_J))return PB11_STATUS_INVALID_ARGUMENT;
 if(cutoff_J<minimum || cutoff_J>maximum)return PB11_STATUS_OUT_OF_RANGE;
 return alpha_impl(l,A,q,cosine,out,cutoff_J/mev,pruned);
}
