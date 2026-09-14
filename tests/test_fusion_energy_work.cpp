#include "fusion_energy_work.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
void check(bool x,const char*m){if(!x)throw std::runtime_error(m);}
void near(double a,double b){check(std::abs(a-b)<=1e-12*std::max(std::abs(a),std::abs(b)),"analytic value mismatch");}
void analytic(){
 double e[]={1,3},old[]={4},out[1],face[2];fusion_energy_work_ledger_v1 l;
 check(fusion_c_energy_work_trial(1,.5,1,e,old,out,face,&l)==0,"cooling single cell");
 near(out[0],3.2);near(l.lower_number_m3,.8);near(l.lower_energy_J_m3,.8);near(l.work_on_particles_J_m3,-.8);
 check(fusion_c_energy_work_trial(1,.5,-1,e,old,out,face,&l)==0,"heating single cell");
 near(out[0],16./7);near(l.upper_number_m3,12./7);near(l.upper_energy_J_m3,36./7);near(l.work_on_particles_J_m3,12./7);
 check(fusion_c_energy_work_trial(1,.5,0,e,old,out,face,&l)==0&&out[0]==4&&l.work_on_particles_J_m3==0,"zero work identity");
 e[0]=0;check(fusion_c_energy_work_trial(1,1,1,e,old,out,face,&l)==0&&out[0]==4&&l.lower_number_m3==0,"zero-energy boundary reflects naturally");
 // This last coarse one-cell state does not resolve cooling; refinement is mandatory.
}
void convergence(){
 for(double c:{-1.,1.}){
  double prior_l1=1e9,prior_u=1e9;
  for(int n:{80,160,320,640}){
   std::vector<double> e(n+1),old(n),out(n),face(n+1);for(int j=0;j<=n;++j)e[j]=8.*j/n;
   for(int j=0;j<n;++j)old[j]=std::max(0.,std::min(e[j+1],4.)-std::max(e[j],2.));
   const int steps=n;const double dt=.2/steps;long double work=0,lossn=0,losse=0;
   fusion_energy_work_ledger_v1 l{};
   for(int k=0;k<steps;++k){
    check(fusion_c_energy_work_trial(n,dt,c,e.data(),old.data(),out.data(),face.data(),&l)==0,"refinement trial");
    work+=l.work_on_particles_J_m3;lossn+=l.lower_number_m3+l.upper_number_m3;losse+=l.lower_energy_J_m3+l.upper_energy_J_m3;
    check(std::all_of(out.begin(),out.end(),[](double x){return std::isfinite(x)&&x>=0;}),"positive finite spectrum");old=out;
   }
   const double s=std::exp(-c*.2);double l1=0,u=0,number=0;
   for(int j=0;j<n;++j){double exact=std::max(0.,std::min(e[j+1],4*s)-std::max(e[j],2*s))/s;l1+=std::abs(out[j]-exact);u+=(e[j]+e[j+1])/2*out[j];number+=out[j];}
   const double uerr=std::abs(u/(6*s)-1);
   check(l1<prior_l1&&uerr<prior_u,"energy/grid time refinement fails");prior_l1=l1;prior_u=uerr;
   check(std::abs(number+lossn-2)<2e-12,"cumulative number ledger");
   check(std::abs(u+losse-6-work)<2e-11,"cumulative work ledger");
   if(n==640)check(uerr<.003,"resolved continuum energy gate");
   std::cout<<"refinement C="<<c<<" cells="<<n<<" steps="<<steps<<" spectrum_L1="<<l1<<" energy_relative="<<uerr<<'\n';
  }
 }
}
void nonuniform_si(){
 double e[]={0,1e-26,1e-18,5e-17,1e-16,1e-15},old[]={1,2,3,4,5},out[5],face[6];
 fusion_energy_work_ledger_v1 l;
 for(double c:{-1.,1.}){
  check(fusion_c_energy_work_trial(5,.2,c,e,old,out,face,&l)==0,"nonuniform SI grid");
  long double ni=0,nf=0,ui=0,uf=0;
  for(int j=0;j<5;++j){ni+=old[j];nf+=out[j];ui+=(static_cast<long double>(e[j])+e[j+1])/2*old[j];uf+=(static_cast<long double>(e[j])+e[j+1])/2*out[j];}
  check(std::abs(nf+l.lower_number_m3+l.upper_number_m3-ni)<1e-12L*ni,"independent SI number balance");
  check(std::abs(uf+l.lower_energy_J_m3+l.upper_energy_J_m3-ui-l.work_on_particles_J_m3)<1e-12L*ui,"independent SI energy balance");
  check(c*l.work_on_particles_J_m3<=0,"signed work on nonuniform grid");
 }
 double scaled[]={1e-16,3e-16},population[]={4};
 check(fusion_c_energy_work_trial(1,.5,1,scaled,population,out,face,&l)==0,"energy scale covariance");near(out[0],3.2);near(l.work_on_particles_J_m3,-.8e-16);
}
void time_refinement(){
 double e[]={1,3},out[1],face[2];fusion_energy_work_ledger_v1 l;
 double previous=1e9;
 for(int steps:{8,16,32,64}){
  double old[]={4};for(int k=0;k<steps;++k){check(fusion_c_energy_work_trial(1,.5/steps,1,e,old,out,face,&l)==0,"temporal trial");old[0]=out[0];}
  double error=std::abs(out[0]-4*std::exp(-.25));check(error<previous,"fixed-grid temporal refinement");previous=error;
 }
 // A wholly subnormal trace follows measured IEEE conversion, not a floor.
 double trace[]={std::numeric_limits<double>::denorm_min()};
 check(fusion_c_energy_work_trial(1,100,1,e,trace,out,face,&l)==0,"subnormal trace work");
}
void invalid(){
 double e[]={0,1,2},old[]={1,2},out[]={9,9},face[]={9,9,9};fusion_energy_work_ledger_v1 l{};
 check(fusion_c_energy_work_trial(2,.1,1,e,old,out,face,&l)==0,"valid setup");
 const double saved0=out[0],saved1=out[1];check(fusion_c_energy_work_trial(2,.1,1,e,old,out,face,&l)==0&&out[0]==saved0&&out[1]==saved1&&old[0]==1&&old[1]==2,"repeatability and input preservation");
 old[1]=-1;check(fusion_c_energy_work_trial(2,.1,1,e,old,out,face,&l)!=0,"negative population accepted");
 check(out[0]==0&&out[1]==0&&face[0]==0&&face[1]==0&&face[2]==0&&l.initial_number_m3==0&&l.work_on_particles_J_m3==0,"failure not cleared");
 old[1]=2;e[2]=1;check(fusion_c_energy_work_trial(2,.1,1,e,old,out,face,&l)!=0,"bad edge accepted");e[2]=2;
 check(fusion_c_energy_work_trial(2,.1,std::numeric_limits<double>::quiet_NaN(),e,old,out,face,&l)!=0,"NaN coefficient accepted");
 check(fusion_c_energy_work_trial(0,1,0,nullptr,nullptr,nullptr,nullptr,nullptr)==PB11_STATUS_INVALID_ARGUMENT,"bad dimensions");
 check(fusion_c_energy_work_trial(1,1,0,nullptr,nullptr,nullptr,nullptr,nullptr)==PB11_STATUS_NULL_OUTPUT,"null outputs");
 old[0]=std::numeric_limits<double>::max();old[1]=old[0];check(fusion_c_energy_work_trial(2,1,0,e,old,out,face,&l)!=0,"ledger overflow accepted");
}
int main(){try{analytic();convergence();nonuniform_si();time_refinement();invalid();std::cout<<"PASS energy work analytic, refinement, signed ledgers and failures\n";return 0;}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
