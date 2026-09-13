#include "fusion_target_burn.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace {
void check(bool x,const char*s){if(!x)throw std::runtime_error(s);}
void close(double a,double b,double tol,const char*s){check(std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)}),s);}
void algebra(){
 double e[]={0,2},N[]={1},K[]={1},M[]={.25},trial[1],loss[1];fusion_target_burn_v1 o{};
 check(fusion_c_target_burn_trial(1,1,e,N,1,1,K,M,trial,loss,&o)==0,"one-bin trial");
 const double exact=(std::sqrt(5.)-1)/2;
 close(trial[0],exact,2e-14,"independent quadratic root");close(o.final_target_number_m3,exact,2e-14,"shared target root");
 close(loss[0],1-exact,2e-14,"particle loss");close(o.final_target_energy_J_m3,1-.25*loss[0],2e-14,"conditioned target energy");
 check(N[0]==1,"old state unchanged");
 double trial2[1],loss2[1];fusion_target_burn_v1 again{};
 check(fusion_c_target_burn_trial(1,1,e,N,1,1,K,M,trial2,loss2,&again)==0 && trial2[0]==trial[0] && loss2[0]==loss[0],"rejected trial has no accumulated state");
 check(fusion_c_target_burn_trial(1,1,e,N,1,.01,K,M,trial,loss,&o)==PB11_STATUS_NUMERICAL_FAILURE && trial[0]==0 && loss[0]==0 && o.reactions_m3==0,"energy overdraw rejects and clears");
 K[0]=0;M[0]=1;
 check(fusion_c_target_burn_trial(1,1,e,N,1,1,K,M,trial,loss,&o)==PB11_STATUS_INVALID_ARGUMENT,"zero rate cannot carry energy");
}
double integrate(int steps){
 double e[]={0,2},N[]={2},K[]={1},M[]={0},trial[1],loss[1],B=1;fusion_target_burn_v1 o{};
 for(int i=0;i<steps;++i){check(fusion_c_target_burn_trial(1,1./steps,e,N,B,0,K,M,trial,loss,&o)==0,"evolving trial");N[0]=trial[0];B=o.final_target_number_m3;}
 close(N[0]-B,1,1e-12,"two-pool invariant over accepted steps");return B;
}
void convergence(){
 const double exact=std::exp(-1.)/(2-std::exp(-1.));
 double a=std::abs(integrate(40)-exact),b=std::abs(integrate(80)-exact),c=std::abs(integrate(160)-exact);
 check(a/b>1.8 && b/c>1.8 && c<.003,"first-order convergence to independent unequal-pool ODE");
 std::cout<<"target ODE errors "<<a<<' '<<b<<' '<<c<<'\n';
}
void multibin(){
 double e[]={0,2,4,8},N[]={10,20,30},K[]={1e10,2e10,0},M[]={2e10,6e10,0},trial[3],loss[3];fusion_target_burn_v1 o{};
 check(fusion_c_target_burn_trial(3,1,e,N,1,10,K,M,trial,loss,&o)==0,"stiff shared target");
 check(o.final_target_number_m3>=0 && o.reactions_m3<=1+1e-12 && trial[2]==30 && loss[2]==0,"no overburn and inert cell unchanged");
 double R=0,removedE=0;for(int i=0;i<3;++i){check(trial[i]>=0 && loss[i]>=0,"positive bins");R+=loss[i];removedE+=loss[i]*(e[i]+e[i+1])/2;}
 close(R,o.reactions_m3,1e-12,"per-bin event ledger");close(removedE,o.removed_fast_energy_J_m3,1e-12,"fast energy debit");
 close(2*loss[0]+3*loss[1],o.removed_target_energy_J_m3,1e-12,"energy-selected target debit");
 close(o.initial_fast_number_m3-o.final_fast_number_m3,R,1e-12,"fast count closure");
 close(o.initial_target_number_m3-o.final_target_number_m3,R,1e-12,"target count closure");
 // Large target reservoir and extreme stiffness exercises the solver's
 // positive lower bound; no artificial fixed relative target cutoff.
 K[0]=1e200;K[1]=1e200;M[0]=0;M[1]=0;
 check(fusion_c_target_burn_trial(3,1,e,N,100,0,K,M,trial,loss,&o)==0,"extremely stiff surplus target");
 close(o.final_target_number_m3,70,1e-12,"all reactive fast bins depleted asymptote");
}
void physical(){
 constexpr double kev=1.602176634e-16,mp=1.67262192595e-27,mb=11*1.66053906892e-27;
 double e[]={400*kev,600*kev},N[]={1e15},trial[1],loss[1];fusion_beam_window_v1 w[1];fusion_target_burn_v1 o{};
 check(fusion_c_beam_target_burn_window_trial(0,1,1.,mp,mb,e,N,1e20,0,trial,loss,w,&o)==0,"physical cold-target window trial");
 check(w[0].resolved_reactivity_m3_s>0 && w[0].domain_incomplete==0 && o.reactions_m3>0,"cold target known-window burning");
 check(o.removed_target_energy_J_m3==0 && o.final_target_energy_J_m3==0,"cold target carries zero initial kinetic energy");
 double residual=(o.removed_fast_energy_J_m3/o.reactions_m3)/(500*kev)-1;
 check(std::abs(residual)<1e-13,"projectile kinetic energy handed to product source");
 check(fusion_c_beam_target_burn_window_trial(0,1,.01,mp,mb,e,N,1e20,10*kev,trial,loss,w,&o)==0,"warm-target window trial");
 check(w[0].domain_incomplete==1,"unresolved nuclear window remains explicit");
 double selected=w[0].target_energy_reactivity_J_m3_s/w[0].resolved_reactivity_m3_s;
 check(std::abs(o.removed_target_energy_J_m3/o.reactions_m3/selected-1)<1e-12,"reaction-conditioned target energy composition");
}
}
int main(){try{algebra();convergence();multibin();physical();std::cout<<"PASS: target burn, conservation, rollback and physical-window composition\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
