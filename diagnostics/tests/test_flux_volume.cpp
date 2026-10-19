#include "fusion_flux_volume_diagnostic.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
const double pi=std::acos(-1.);
void require(bool ok,const char* what){if(!ok)throw std::runtime_error(what);}
void near(double actual,double expected,double tol=3e-12){
 require(std::abs(actual-expected)<=tol*std::max(std::abs(expected),1e-25),"volume mismatch");
}
struct Snapshot {
 int m; double shift, elong;
 std::vector<double> x{0,.2,.5,1},a,d,F{1.6,1.6,1.6,1.6},C{0,.016,.04,.08};
 Snapshot(int modes=1,double c=0,double b=.6):m(modes),shift(c),elong(b),a(4*(2*m+1)),d(a.size()){
  for(int i=0;i<4;++i){double q=x[i];int k=i*(2*m+1);
   a[k]=2+c*q*q;d[k]=2*c*q;
   for(int h=1;h<=m;++h){double ar=h==1?.4:.012/h,az=h==1?b:.009/h;
    a[k+h]=ar*std::pow(q,h);a[k+m+h]=az*std::pow(q,h);
    d[k+h]=ar*h*std::pow(q,h-1);d[k+m+h]=az*h*std::pow(q,h-1);
   }
  }
 }
 int call(double q,fusion_flux_volume_diagnostic_v1* out,int core=1)const{
  return fusion_c_flux_core_enclosed_volume_diagnostic_v1(4,m,x.data(),a.data(),d.data(),F.data(),C.data(),core,q,out);
 }
 double analytic(double q)const{return 2*pi*pi*(2+shift*q*q)*.4*elong*q*q;}
 // Independent dense midpoint line integral of the explicitly prescribed shape;
 // this never calls either diagnostic or field evaluator and keeps the full R^2.
 double dense(double q)const{
  constexpr int count=65536;long double sum=0;
  for(int j=0;j<count;++j){long double t=2*static_cast<long double>(pi)*(j+.5L)/count;
   long double R=2+shift*q*q,Zt=0;
   for(int h=1;h<=m;++h){double ar=h==1?.4:.012/h,az=h==1?elong:.009/h;
    R+=ar*std::pow(q,h)*std::cos(h*t);Zt+=h*az*std::pow(q,h)*std::cos(h*t);
   }
   sum+=R*R*Zt;
  }
  return static_cast<double>(2*static_cast<long double>(pi)*pi*sum/count);
 }
};
bool cleared(const fusion_flux_volume_diagnostic_v1& v){
 return v.version==0&&v.angular_nodes==0&&v.boundary_orientation==0&&v.reserved==0&&
 v.signed_volume_m3==0&&v.absolute_volume_m3==0&&v.min_abs_boundary_jacobian_m3==0;
}
void poison(fusion_flux_volume_diagnostic_v1&v){v={7,7,7,7,7,7,7};}
}
int main(){try{
 fusion_flux_volume_diagnostic_v1 v{};
 for(double elong:{.4,.6,-.6})for(double shift:{0.,.025}){
  Snapshot s(1,shift,elong);
  for(int core:{1,2}){
   require(s.call(0,&v,core)==PB11_STATUS_OK,"axis status");
   require(v.version==1&&v.signed_volume_m3==0&&v.absolute_volume_m3==0&&v.boundary_orientation==0,"axis exact zero");
   double previous=0;
   for(double q:{1e-6,.03,.199999,.2,.200001,.499999,.5,.500001,.8,1.}){
    require(s.call(q,&v,core)==PB11_STATUS_OK,"analytic status");
    near(v.signed_volume_m3,s.analytic(q));near(v.absolute_volume_m3,std::abs(s.analytic(q)));
    require(v.version==1&&v.reserved==0&&v.angular_nodes>3&&v.min_abs_boundary_jacobian_m3>0,"metadata");
    require(v.boundary_orientation==(elong>0?1:-1),"orientation");
    double shell=v.absolute_volume_m3-previous;
    require(shell>0,"positive shell difference");previous=v.absolute_volume_m3;
   }
   for(double q:{s.x[core]}){
    fusion_flux_volume_diagnostic_v1 left{},right{};
    require(s.call(std::nextafter(q,0.),&left,core)==0&&s.call(std::nextafter(q,1.),&right,core)==0,"matching status");
    near(left.signed_volume_m3,right.signed_volume_m3);
   }
   fusion_flux_volume_diagnostic_v1 inner{},outer{};
   require(s.call(.15,&inner,core)==0&&s.call(.85,&outer,core)==0,"shell status");
   near(outer.signed_volume_m3-inner.signed_volume_m3,s.analytic(.85)-s.analytic(.15));
  }
 }
 for(int m:{2,3}){Snapshot s(m,.015);
  for(int core:{1,2})for(double q:{.07,.2,.37,.5,.9,1.}){
   require(s.call(q,&v,core)==0,"multiharmonic status");near(v.signed_volume_m3,s.dense(q),2e-12);
   require(v.angular_nodes>3*m&&(v.angular_nodes&(v.angular_nodes-1))==0,"angular rule");
  }
 }
 Snapshot maximum(64);
 for(int i=0;i<4;++i)for(int h=2;h<=64;++h){
  const int k=i*129;
  maximum.a[k+h]=maximum.a[k+64+h]=0;
  maximum.d[k+h]=maximum.d[k+64+h]=0;
 }
 for(double q:{.37,1.}){
  require(maximum.call(q,&v)==PB11_STATUS_OK,"maximum harmonics status");
  near(v.signed_volume_m3,maximum.analytic(q));
  require(v.angular_nodes==256,"maximum harmonics angular rule");
 }
 Snapshot s;
 auto invalid=[&](int n,int m,const double*x,const double*a,const double*d,const double*F,const double*C,int core,double q,int expected){
  poison(v);int rc=fusion_c_flux_core_enclosed_volume_diagnostic_v1(n,m,x,a,d,F,C,core,q,&v);
  require(rc==expected,"invalid-input status");require(cleared(v),"error output cleared");
 };
 auto bad=[&](double q,int expected){invalid(4,1,s.x.data(),s.a.data(),s.d.data(),s.F.data(),s.C.data(),1,q,expected);};
 require(s.call(.3,nullptr)==PB11_STATUS_NULL_OUTPUT,"null output");
 for(double q:{-.01,1.01})bad(q,PB11_STATUS_OUT_OF_RANGE);
 bad(std::numeric_limits<double>::quiet_NaN(),PB11_STATUS_INVALID_ARGUMENT);
 for(int which=0;which<5;++which)
  invalid(4,1,which==0?nullptr:s.x.data(),which==1?nullptr:s.a.data(),which==2?nullptr:s.d.data(),which==3?nullptr:s.F.data(),which==4?nullptr:s.C.data(),1,.3,PB11_STATUS_INVALID_ARGUMENT);
 for(int n:{-1,0,2,10001})invalid(n,1,s.x.data(),s.a.data(),s.d.data(),s.F.data(),s.C.data(),1,.3,PB11_STATUS_OUT_OF_RANGE);
 for(int m:{0,65})invalid(4,m,s.x.data(),s.a.data(),s.d.data(),s.F.data(),s.C.data(),1,.3,PB11_STATUS_OUT_OF_RANGE);
 for(int core:{0,3})invalid(4,1,s.x.data(),s.a.data(),s.d.data(),s.F.data(),s.C.data(),core,.3,PB11_STATUS_OUT_OF_RANGE);
 for(auto* values:{&s.x,&s.a,&s.d,&s.F,&s.C}){
  double saved=(*values)[1];(*values)[1]=std::numeric_limits<double>::infinity();
  bad(.3,PB11_STATUS_INVALID_ARGUMENT);(*values)[1]=saved;
 }
 Snapshot folded(1,1.);poison(v);
 require(folded.call(1,&v)==PB11_STATUS_NUMERICAL_FAILURE&&cleared(v),"boundary orientation failure");
 std::cout<<"flux volume analytic, dense Fourier, continuity, shell and invalid-input checks passed\n";
 return EXIT_SUCCESS;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return EXIT_FAILURE;}
 catch(...){std::cerr<<"unexpected exception escaped diagnostic tests\n";return EXIT_FAILURE;}
}
