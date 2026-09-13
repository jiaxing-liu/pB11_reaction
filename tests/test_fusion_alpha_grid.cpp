#include "fusion_alpha_spectrum.h"
#include "fusion_alpha_amplitudes.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
namespace {
constexpr double mev=1.602176634e-13;
void require(bool v,const char*why){if(!v)throw std::runtime_error(why);}
double amplitude_norm(const fusion_alpha_amplitudes_v1&a){double s=0;for(int i=0;i<5;++i)s+=a.sym_real[i]*a.sym_real[i]+a.sym_imag[i]*a.sym_imag[i];return s;}
struct Result{std::vector<double> birth;fusion_alpha_spectrum_v1 ledger;};
std::vector<double> edges(int n,double lo=0,double hi=8){std::vector<double> r(n+1);for(int i=0;i<=n;++i)r[i]=(lo+(hi-lo)*i/n)*mev;return r;}
Result evaluate(int mode,int n,double k=.76,double phase=.67*2*3.141592653589793,double cutoff=.001,const std::vector<double>&grid=edges(32)){
 Result r;r.birth.resize(grid.size()-1);
 int status=fusion_c_alpha_spectrum_grid(mode,9.3*mev,cutoff*mev,k,phase,n,n,int(r.birth.size()),grid.data(),r.birth.data(),&r.ledger);
 if(status)std::cerr<<"mode="<<mode<<" n="<<n<<" status="<<status<<'\n';
 require(status==0,"grid status");
 const auto &l=r.ledger;
 require(std::abs(l.mapped_number+l.below_number+l.above_number-3)<4e-12,"three alpha closure including spill");
 require(std::abs((l.mapped_energy_J+l.below_energy_J+l.above_energy_J)/mev-9.3)<2e-11,"CM energy closure including spill");
 double number=0,energy=0;for(std::size_t i=0;i<r.birth.size();++i){require(r.birth[i]>=0 && std::isfinite(r.birth[i]),"finite nonnegative births");number+=r.birth[i];energy+=r.birth[i]*(grid[i]+grid[i+1])/2;}
 require(std::abs(number-l.mapped_number)<4e-12,"mapped count diagnostic");
 require(std::abs((energy-l.mapped_energy_J)/mev)<2e-11,"mapped energy diagnostic");
 require(l.normalization_J2>0 && l.quadrature_events==n*n,"normalization and quadrature count");return r;
}
double distance(const Result&a,const Result&b){double s=0;for(std::size_t i=0;i<a.birth.size();++i)s+=std::abs(a.birth[i]-b.birth[i]);return s/3;}
void cutoff(){
 fusion_alpha_amplitudes_v1 strict{},cut{};int pruned=-1;
 require(fusion_c_alpha_amplitudes(1,9.3*mev,3.129*mev,.3,&strict)==0,"strict event");
 require(fusion_c_alpha_amplitudes_cutoff(1,9.3*mev,3.129*mev,.3,.001*mev,&cut,&pruned)==0 && pruned==0,"ordinary unpruned event");
 for(int i=0;i<5;++i)require(strict.sym_real[i]==cut.sym_real[i] && strict.sym_imag[i]==cut.sym_imag[i],"strict API preserved exactly");
 // Even primary l admits nonzero collinear endpoint amplitudes; odd l
 // vanishes there by the J=2 angular selection rule.
 for(double q:{0.,9.3}){
  require(fusion_c_alpha_amplitudes_cutoff(2,9.3*mev,q*mev,.3,.001*mev,&cut,&pruned)==0 && pruned==1,"only endpoint permutation is pruned");
  require(amplitude_norm(cut)>0 && cut.phase_space_J==0,"endpoint retains other amplitudes");
 }
 require(fusion_c_alpha_amplitudes_cutoff(1,9.3*mev,3.129*mev,.3,.0009*mev,&cut,&pruned)==PB11_STATUS_OUT_OF_RANGE && pruned==0 && cut.phase_space_J==0,"cutoff domain clears output");
}
void grid(){
 auto one=evaluate(1,40),three=evaluate(3,40),mixed_one=evaluate(13,40,1),mixed_three=evaluate(13,40,0);
 require(distance(one,mixed_one)<1e-13 && distance(three,mixed_three)<1e-13,"unit-normalized mixture endpoint shapes");
 auto phase0=evaluate(13,40,.76,0),phasepi=evaluate(13,40,.76,3.141592653589793);
 require(distance(phase0,phasepi)>1e-3,"coherent phase changes spectrum");
 auto spill=evaluate(2,40,.76,0,.001,edges(20,2,4));
 require(spill.ledger.below_number>0 && spill.ledger.above_number>0,"spill retained on both sides");
 auto strictcut=evaluate(13,64),largecut=evaluate(13,64,.76,.67*2*3.141592653589793,.01);
 require(distance(strictcut,largecut)<1e-9,"sampled cutoff sensitivity");
 auto fine=evaluate(13,128),coarse=evaluate(13,64);
 double difference=distance(fine,coarse);
 std::cout<<"64/128 per-event histogram normalized L1 difference "<<difference<<'\n';
 require(difference<.02,"representative grid refinement");
 std::vector<double> e=edges(2),birth(2,1.);fusion_alpha_spectrum_v1 l{};l.mapped_number=1;
 int status=fusion_c_alpha_spectrum_grid(13,9.3*mev,.001*mev,.76,0,3,4,2,e.data(),birth.data(),&l);
 require(status==PB11_STATUS_INVALID_ARGUMENT && birth[0]==0 && birth[1]==0 && l.mapped_number==0,"invalid quadrature clears outputs");
 e[1]=std::numeric_limits<double>::quiet_NaN();birth[0]=1;
 status=fusion_c_alpha_spectrum_grid(2,9.3*mev,.001*mev,.76,0,4,4,2,e.data(),birth.data(),&l);
 require(status==PB11_STATUS_INVALID_ARGUMENT && birth[0]==0,"invalid grid rejected");
}
}
int main(){try{cutoff();grid();std::cout<<"PASS: explicit cutoff, normalized spectrum, coherent mixture, spill and refinement\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
