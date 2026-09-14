#include "fusion_target_network.h"
#include "fusion_beam_birth.h"
#include "fusion_nuclear_data.h"
#include <array>
#include <vector>
#include <cmath>
#include <cstdio>
int main(){constexpr double keV=1.602176634e-16,MeV=1.602176634e-13;constexpr int n=64;
 std::vector<double> grid(n+1);for(int j=0;j<=n;++j)grid[j]=25*MeV*j/n;
 fusion_beam_birth_options_v1 options{2.5*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,16,8,8};
 double E[2]={100*keV,80*keV},F[2]={1e18,2e18},B[2]={1e19,5e18},U[2]={1.5e19*10*keV,1.5*5e18*10*keV},nf[2],nb[2],nu[2],loss[4];
 int ch[4]={1,2,3,3},slot[4]={0,0,0,1},fi[4]={0,0,0,1},ti[4]={0,0,1,0};
 std::array<fusion_beam_birth_v1,4> source{};std::array<fusion_target_network_edge_v1,4> edges{};std::array<std::vector<double>,4> spectra;
 for(int e=0;e<4;++e){spectra[e].resize(7*n);int st=fusion_c_beam_birth_grid(ch[e],slot[e],E[fi[e]],10*keV,&options,n,grid.data(),spectra[e].data(),&source[e]);if(st)return st;auto&s=source[e].spectrum;if(!(s.reactivity_m3_s>0))return 1;edges[e]={fi[e],ti[e],s.reactivity_m3_s,s.reactant_energy_moment_J_m3_s[1-slot[e]]};}
 fusion_target_network_v1 result{};int st=fusion_c_target_network_trial(2,2,4,1,E,F,B,U,edges.data(),nf,nb,nu,loss,&result);if(st)return st;
 long double number[7]{},energy=0,Q=0,expected[7]{};
 for(int e=0;e<4;++e){const auto&s=source[e].spectrum;long double scale=static_cast<long double>(loss[e])/s.reactivity_m3_s;fusion_nuclear_channel_v1 reaction{};if(fusion_c_nuclear_channel(ch[e],&reaction))return 1;Q+=static_cast<long double>(loss[e])*reaction.q_J;for(int j=0;j<reaction.product_count;++j)expected[reaction.product_ids[j]]+=loss[e];
  for(int id=0;id<7;++id){number[id]+=scale*(s.below_number_m3_s[id]+s.above_number_m3_s[id]);energy+=scale*(s.below_energy_J_m3_s[id]+s.above_energy_J_m3_s[id]);for(int j=0;j<n;++j){number[id]+=scale*spectra[e][id*n+j];energy+=scale*spectra[e][id*n+j]*(grid[j]+grid[j+1])/2;}}
 }
 long double maxN=0;for(int id=0;id<7;++id){long double err=expected[id]>0?std::abs(number[id]/expected[id]-1):std::abs(number[id]);if(!std::isfinite(err))return 1;maxN=std::max(maxN,err);}
 long double budget=Q+result.removed_fast_energy_J_m3+result.removed_target_energy_J_m3;long double er=energy/budget-1;
 // DD branch ratio must equal K1/K2: both channels use the same final pools.
 long double branch=(static_cast<long double>(loss[0])/loss[1])/(edges[0].reactivity_m3_s/edges[1].reactivity_m3_s)-1;
 std::printf("events %.17g max_species_relative %.12Lg energy_relative %.12Lg DD_branch_relative %.12Lg iterations %d\n",result.reactions_m3,maxN,er,branch,result.iterations);
 return !std::isfinite(er)||!std::isfinite(branch)||maxN>3e-10L||std::abs(er)>3e-10L||std::abs(branch)>3e-12L;
}
