#include "fusion_radial_transport.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>
void check(bool x,const char*m){if(!x)throw std::runtime_error(m);}
void near(double a,double b){check(std::abs(a-b)<=2e-13*std::max(std::abs(a),std::abs(b)),"analytic mismatch");}
struct Run {
 int z,c;std::vector<double> vo,vn,a,k,n,bc,out,face;std::vector<fusion_radial_ledger_v1> ledger;
 Run(int zones,int comps):z(zones),c(comps),vo(z,1),vn(z,1),a(z+1),k(c*(z+1)),n(c*z),bc(2*c),out(c*z),face(c*(z+1)),ledger(c){}
 int call(double dt){return fusion_c_radial_transport_trial(z,c,dt,vo.data(),vn.data(),a.data(),k.data(),n.data(),bc.data(),out.data(),face.data(),ledger.data());}
 void balance(){for(int p=0;p<c;++p){long double old=0,next=0;for(int j=0;j<z;++j){check(out[p*z+j]>=0,"negative output");old+=static_cast<long double>(vo[j])*n[p*z+j];next+=static_cast<long double>(vn[j])*out[p*z+j];}long double r=next-old-face[p*(z+1)]+face[p*(z+1)+z];check(std::abs(r)<1e-12L*(std::abs(old)+std::abs(next)+std::abs(face[p*(z+1)])+std::abs(face[p*(z+1)+z])),"independent extensive balance");}}
 void clear(){check(std::all_of(out.begin(),out.end(),[](double x){return x==0;}),"density not cleared");check(std::all_of(face.begin(),face.end(),[](double x){return x==0;}),"faces not cleared");for(auto l:ledger)check(l.initial_number==0&&l.final_number==0&&l.inner_inward_number==0&&l.outer_outward_number==0&&l.balance_error==0,"ledger not cleared");}
};
void test_geometry(){
 Run r(2,2);r.vo={1,2};r.vn={1.1,1.9};r.a={0,-1,0};r.n={3,3,7,7};
 const auto input=r.n;check(r.call(.1)==0,"moving-grid uniform state");
 for(int i=0;i<4;++i)near(r.out[i],r.n[i]);near(r.face[1],-.3);near(r.face[4],-.7);r.balance();
 const auto oldout=r.out,oldface=r.face;check(r.call(.1)==0&&r.out==oldout&&r.face==oldface&&r.n==input,"trial repeatability/input preservation");
 r.a={0,0,0};r.vn={.5,1};check(r.call(.1)==0,"compression dilution");for(int i=0;i<4;++i)near(r.out[i],2*r.n[i]);r.balance();
}
void test_diffusion(){
 Run r(2,2);r.n={10,0,0,20};r.k={0,2,0,0,1,0};check(r.call(.1)==0,"two species diffusion");
 near(r.out[0],60./7);near(r.out[1],10./7);near(r.out[2],5./3);near(r.out[3],55./3);near(r.face[1],10./7);near(r.face[4],-5./3);r.balance();
 // Independent spectral solution for a symmetric two-cell exchange.
 double prior=1e99;
 for(int steps:{8,16,32,64}){Run t(2,1);t.n={10,0};t.k={0,2,0};for(int i=0;i<steps;++i){check(t.call(.2/steps)==0,"diffusion refinement");t.n=t.out;}
 double exact=5+5*std::exp(-.8),error=std::abs(t.out[0]-exact);check(error<prior,"time refinement did not improve");prior=error;}
}
void test_boundaries(){
 Run r(1,1);r.n={10};r.bc={3,0};r.a={2,1};check(r.call(.5)==0,"open advection");near(r.out[0],26./3);near(r.face[0],3);near(r.face[1],13./3);r.balance();
 r.bc={0,4};r.a={-2,-1};check(r.call(.5)==0,"reversed advection");near(r.out[0],6);near(r.face[0],-6);near(r.face[1],-2);r.balance();
 r.n={0};r.a={0,0};r.k={1,1};r.bc={1,3};check(r.call(1)==0,"diffusive reservoirs");near(r.out[0],4./3);r.balance();
}
void test_validation(){
 Run r(2,2);r.n={1,2,3,4};check(r.call(1)==0,"initial valid call");r.k.back()=-1;check(r.call(1)!=0,"last component invalid accepted");r.clear();r.k.back()=0;
 r.vo[0]=0;check(r.call(1)!=0,"zero volume accepted");r.clear();r.vo[0]=1;
 r.a[1]=std::numeric_limits<double>::quiet_NaN();check(r.call(1)!=0,"NaN accepted");r.clear();r.a[1]=0;
 check(r.call(0)!=0,"zero dt accepted");r.clear();
 r.n[0]=-1;check(r.call(1)!=0,"negative density accepted");r.clear();r.n[0]=1;
 r.vo[0]=std::numeric_limits<double>::max();r.n[0]=std::numeric_limits<double>::max();check(r.call(1)!=0,"overflow accepted");r.clear();
 check(fusion_c_radial_transport_trial(0,1,1,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr)==PB11_STATUS_INVALID_ARGUMENT,"invalid dimensions");
 check(fusion_c_radial_transport_trial(1,1,1,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr)==PB11_STATUS_NULL_OUTPUT,"null outputs");
 Run tiny(1,1);tiny.vo={std::numeric_limits<double>::denorm_min()};tiny.vn={1};tiny.n={.125};check(tiny.call(1)==0&&tiny.out[0]==0,"measured subnormal rounding");
}
int main(){try{test_geometry();test_diffusion();test_boundaries();test_validation();std::cout<<"PASS radial transport: geometry, analytic fluxes, positivity, balance and failure atomicity\n";return 0;}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}}
