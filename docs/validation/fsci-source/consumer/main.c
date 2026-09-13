#include <fusion_alpha_spectrum.h>
#include <fusion_alpha_amplitudes.h>
#include <fusion_nuclear_coulomb.h>
#include <math.h>
int main(void){const double mev=1.602176634e-13;double edges[21],birth[20];
 fusion_nuclear_coulomb_v1 coulomb={0};fusion_alpha_amplitudes_v1 amplitude={0};fusion_alpha_spectrum_v1 out={0};int pruned=0;
 for(int i=0;i<=20;++i)edges[i]=i*.5*mev;
 if(fusion_c_nuclear_coulomb_radius16(1,mev,&coulomb))return 1;
 if(fusion_c_alpha_amplitudes_fsci_cutoff(2,1,8.84*mev,3.129*mev,0,.001*mev,&amplitude,&pruned))return 2;
 if(fusion_c_alpha_spectrum_model_grid(2,1,8.84*mev,.001*mev,.76,4.2,16,16,20,edges,birth,&out))return 3;
 if(fabs(out.mapped_number+out.below_number+out.above_number-3)>1e-10)return 4;
 return 0;}
