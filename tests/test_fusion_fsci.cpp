#include "fusion_nuclear_coulomb.h"
#include "fusion_alpha_amplitudes.h"
#include "fusion_alpha_spectrum.h"
#include <cmath>
#include <algorithm>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>
constexpr double mev=1.602176634e-13;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void ok(int value){require(value==0,"API failed");}
int main(){try{
 // Independent mpmath40-digit fixed-radius reference points; see fsci-r16-validation.json.
 const double coulomb_refs[][6]={
  {0,0.001,-388.74056520001488,-4.7984037314727583,-0.06082775669595343,0.99814827756968949},
  {0,0.002,-268.99230713723642,-4.7948152295793038,0.85073025168908545,0.52560254837759812},
  {0,0.01,-109.21330411324409,-4.7659929843970064,0.99968210736603635,-0.025212778748100004},
  {0,0.10000000000000001,-21.228663630612942,-4.4261697578890056,-0.017174018288112297,0.99985251567210631},
  {0,1,1.1011347823113724,-0.53462133006966983,0.9696633721662975,0.24444415452017773},
  {0,4,2.3285381368140365,-0.064038418089615926,0.87959087318260465,-0.47573090693569942},
  {0,12,2.9521505212498349,-0.018676241728666641,-0.92217947255396504,0.38676222721472031},
  {1,0.001,-390.27607584482178,-5.1507574351917942,-0.99954582512891699,-0.030135418818267635},
  {1,0.002,-270.52630837272073,-5.1474067965668517,-0.48816388482713852,0.87275198169381218},
  {1,0.01,-110.73523809589209,-5.1205083027049429,0.12180061292162628,0.99255458826802878},
  {1,0.10000000000000001,-22.615261034403453,-4.8053429359691266,-0.96077540678244977,0.27732763605890082},
  {1,1,0.92932326573442647,-0.79392918165017823,0.18000539277265556,0.98366562335621044},
  {1,4,2.3095520037202948,-0.085389264483661048,0.6283693626969642,0.77791512649119443},
  {1,12,2.9466898918545268,-0.024344873201095935,-0.47917832213082145,-0.87771757165952347},
  {2,0.001,-392.52304909296709,-5.6446525450154521,-0.015929253764621044,-0.99987312138816009},
  {2,0.002,-272.77060772712434,-5.6415907435155237,-0.90266316806790947,-0.43034777217223397},
  {2,0.01,-112.95828491653405,-5.6170248779408194,-0.96460028086265948,0.26371632137522044},
  {2,0.10000000000000001,-24.613574250146126,-5.3312008038807956,-0.65410302647402152,-0.75640546716528001},
  {2,1,0.61287851837238638,-1.2988175112300147,-0.6723961742817508,0.74019145145784093},
  {2,4,2.2797460973914547,-0.12035651185710515,-0.57602492402512673,0.81743212984433578},
  {2,12,2.9383873460014835,-0.033079189600099097,0.79368041485527541,-0.60833493987700382},
  {3,0.001,-164.19685731517674,-2.7425029292824816,-0.76381442127935972,-0.64543592233906277},
  {3,0.002,-112.34932915919067,-2.7376524446359034,-0.93174511099856838,0.36311299636651068},
  {3,0.01,-43.203238582451,-2.6984429222218878,0.35139583508388783,-0.93622698480961186},
  {3,0.10000000000000001,-5.6144284469368833,-2.1827472178736236,0.65530696637051622,0.75536268098594272},
  {3,1,1.3854557336768205,-0.12321738984088601,-0.98719378707781813,-0.15952563040137299},
  {3,4,2.2457334217402272,-0.024273751809385017,-0.66349220369507345,-0.74818319657411125},
  {3,12,2.8267083489637734,-0.0076899622675898347,-0.86130981198397616,0.50808011944980447},
  {4,0.001,-167.79307337394252,-3.5290490207013958,0.8280036583695094,0.56072269592616697},
  {4,0.002,-115.93687140499996,-3.5252412771461779,0.86692349243369504,-0.49844122849796957},
  {4,0.01,-46.722082529608819,-3.494579474563289,-0.026500317859774941,0.99964880490766905},
  {4,0.10000000000000001,-8.4059715623955942,-3.1188409055410697,-0.99486266000788082,0.10123382696531708},
  {4,1,1.1945207600614713,-0.35605779351820893,0.87759863937469129,-0.47939610779363923},
  {4,4,2.2116485604327152,-0.060444935657084836,0.86006235692318278,0.51018892795095017},
  {4,12,2.8161392931657256,-0.018473030139908777,0.75887300742580299,-0.65123863414305916}
 };
 for(const auto& row:coulomb_refs){fusion_nuclear_coulomb_v1 c{};
  ok(fusion_c_nuclear_coulomb_radius16(int(row[0]),row[1]*mev,&c));
  require(std::abs(c.log_penetrability-row[2])<1e-7,"radius16 logP reference");
  require(std::abs(c.shift-row[3])<1e-8*std::max(1.,std::abs(row[3])),"radius16 shift reference");
  require(std::hypot(c.phase_real-row[4],c.phase_imag-row[5])<1.1e-7,"radius16 phase reference");
 }
 // Direct mpmath50-digit references use Coulomb F/G at16fm, not generated tables.
 const double refs[][5]={{1,8.84,3.129,0,.97397261197406232562},
 {2,8.84,3.129,.75,.95054821814986487899},
 {3,9.3,3.129,-.3,.98643694401709908972},
 {2,8.84,.1,.8,.97902537242645326321},
 {2,8.84,8.6,-.4,106.73940890691405162},
 {2,8.84,4.42,0,.98409209034897500746}};
 for(const auto& row:refs){fusion_alpha_amplitudes_v1 ordinary{},corrected{},none{};int p=0;
  ok(fusion_c_alpha_amplitudes_cutoff(int(row[0]),row[1]*mev,row[2]*mev,row[3],.001*mev,&ordinary,&p));
  ok(fusion_c_alpha_amplitudes_fsci_cutoff(int(row[0]),1,row[1]*mev,row[2]*mev,row[3],.001*mev,&corrected,&p));
  require(p==0,"unexpected suppression in interior reference");
  ok(fusion_c_alpha_amplitudes_fsci_cutoff(int(row[0]),0,row[1]*mev,row[2]*mev,row[3],.001*mev,&none,&p));
  double scale=0,error=0;
  for(int j=0;j<5;++j){
   std::complex<double> a(ordinary.unsym_real[j],ordinary.unsym_imag[j]);
   std::complex<double> b(corrected.unsym_real[j],corrected.unsym_imag[j]);
   scale+=std::norm(row[4]*a);error+=std::norm(b-row[4]*a);
   require(none.sym_real[j]==ordinary.sym_real[j]&&none.sym_imag[j]==ordinary.sym_imag[j],"NONE compatibility");
  }
  require(scale>0&&std::sqrt(error/scale)<5e-7,"direct finite-radius amplitude reference");
 }
 const int n=40;std::vector<double> edges(n+1),birth(n),old(n);
 for(int i=0;i<=n;++i)edges[i]=10*mev*i/n;
 for(int mode:{1,2,3,13}){fusion_alpha_spectrum_v1 a{},b{};
  ok(fusion_c_alpha_spectrum_model_grid(mode,1,8.84*mev,.001*mev,.76,4.2,32,32,n,edges.data(),birth.data(),&a));
  require(std::abs(a.mapped_number+a.below_number+a.above_number-3)<1e-10,"three-alpha number");
  require(std::abs(a.mapped_energy_J+a.below_energy_J+a.above_energy_J-8.84*mev)<1e-10*8.84*mev,"CM source energy");
  for(auto value:birth)require(std::isfinite(value)&&value>=0,"positive source");
  ok(fusion_c_alpha_spectrum_grid(mode,8.84*mev,.001*mev,.76,4.2,16,16,n,edges.data(),old.data(),&b));
  ok(fusion_c_alpha_spectrum_model_grid(mode,0,8.84*mev,.001*mev,.76,4.2,16,16,n,edges.data(),birth.data(),&a));
  require(old==birth,"NONE grid compatibility");
 }
 fusion_nuclear_coulomb_v1 c{},original{};
 ok(fusion_c_nuclear_coulomb_radius16(1,mev,&c));ok(fusion_c_nuclear_coulomb(1,mev,&original));
 require(std::abs(c.rho/original.rho-16/5.1)<1e-12,"explicit radius");
 require(std::abs(c.log_penetrability-original.log_penetrability)>1e-3,"not merely changed rho");
 require(fusion_c_nuclear_coulomb_radius16(1,0,&c)!=0&&c.rho==0,"Coulomb domain clearing");
 fusion_alpha_amplitudes_v1 a{};int pruned=9;
 require(fusion_c_alpha_amplitudes_fsci_cutoff(2,7,8.84*mev,3*mev,0,.001*mev,&a,&pruned)!=0&&pruned==0&&a.phase_space_J==0,"policy rejection");
 ok(fusion_c_alpha_amplitudes_fsci_cutoff(2,1,8.84*mev,0,0,.001*mev,&a,&pruned));
 require(pruned==3,"additional pair cutoff must apply in every permutation");
 for(auto v:a.sym_real)require(v==0,"suppressed real amplitude");
 for(auto v:a.sym_imag)require(v==0,"suppressed imaginary amplitude");
 fusion_alpha_spectrum_v1 grid{};
 require(fusion_c_alpha_spectrum_model_grid(2,7,8.84*mev,.001*mev,.76,4.2,16,16,n,edges.data(),birth.data(),&grid)!=0,"invalid grid policy");
 for(auto v:birth)require(v==0,"invalid grid output clearing");
 require(fusion_c_alpha_spectrum_model_grid(2,1,8.84*mev,.001*mev,.76,4.2,1025,16,n,edges.data(),birth.data(),&grid)!=0,"quadrature upper bound");
 std::cout<<"All finite-radius FSCI tests passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
