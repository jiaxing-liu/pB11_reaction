#include <fusion_beam_birth.h>
#include <fusion_target_burn.h>
#include <fusion_nuclear_data.h>
#include <array>
#include <vector>
#include <cmath>
#include <cstdio>
int main(){constexpr double k=1.602176634e-16,m=1.602176634e-13;const int n=32;
std::vector<double> edges(n+1);for(int j=0;j<=n;++j)edges[j]=25*m*j/n;
fusion_beam_birth_options_v1 o{2.5*m,40,.09184*m,.001*m,.76,0,.051,1,1,0,0,13,0,16,16,8,8};
std::array<double,3> fast_edges{25*k,75*k,125*k};std::array<double,2> old{1e17,2e17},trial{},loss{},K{},Mt{};
std::array<std::vector<double>,2> grids{std::vector<double>(7*n),std::vector<double>(7*n)};std::array<fusion_beam_birth_v1,2> sources{};
for(int i=0;i<2;++i){double E=(fast_edges[i]+fast_edges[i+1])/2;int s=fusion_c_beam_birth_grid(3,0,E,10*k,&o,n,edges.data(),grids[i].data(),&sources[i]);if(s)return s;K[i]=sources[i].spectrum.reactivity_m3_s;if(!(K[i]>0))return 1;Mt[i]=sources[i].spectrum.reactant_energy_moment_J_m3_s[1];}
fusion_target_burn_v1 burn{};int s=fusion_c_target_burn_trial(2,1.,fast_edges.data(),old.data(),1e19,1.5*1e19*10*k,K.data(),Mt.data(),trial.data(),loss.data(),&burn);if(s)return s;
long double number=0,energy=0;for(int i=0;i<2;++i){long double scale=static_cast<long double>(loss[i])/K[i];for(int id=0;id<7;++id){auto&r=sources[i].spectrum;number+=scale*(r.below_number_m3_s[id]+r.above_number_m3_s[id]);energy+=scale*(r.below_energy_J_m3_s[id]+r.above_energy_J_m3_s[id]);for(int j=0;j<n;++j){number+=scale*grids[i][id*n+j];energy+=scale*grids[i][id*n+j]*(edges[j]+edges[j+1])/2;}}}
fusion_nuclear_channel_v1 reaction{};fusion_c_nuclear_channel(3,&reaction);long double expected=static_cast<long double>(reaction.q_J)*burn.reactions_m3+burn.removed_fast_energy_J_m3+burn.removed_target_energy_J_m3;
long double nr=number/(2*burn.reactions_m3)-1,er=energy/expected-1;
std::printf("events=%.17g target_final=%.17g removed_fast_J_m3=%.17g removed_target_J_m3=%.17g number_relative=%.12Lg energy_relative=%.12Lg\n",burn.reactions_m3,burn.final_target_number_m3,burn.removed_fast_energy_J_m3,burn.removed_target_energy_J_m3,nr,er);
return !std::isfinite(nr)||!std::isfinite(er)||std::abs(nr)>3e-10||std::abs(er)>3e-10||burn.reactions_m3<=0;}
