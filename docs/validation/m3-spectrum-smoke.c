#include "fusion_nuclear_coulomb.h"
#include "fusion_alpha_amplitudes.h"
#include <math.h>
#include <stdio.h>
int main(void){
 const double mev=1.602176634e-13;
 fusion_nuclear_coulomb_v1 c;fusion_alpha_amplitudes_v1 a;
 if(fusion_c_nuclear_coulomb(FUSION_ALPHA_ALPHA_L2,3.129*mev,&c))return 1;
 if(fabs(c.shift+.923134650192041157)>1e-7)return 2;
 if(fusion_c_alpha_amplitudes(2,9.3*mev,3.129*mev,.3,&a))return 3;
 double w=0;for(int i=0;i<5;++i)w+=a.sym_real[i]*a.sym_real[i]+a.sym_imag[i]*a.sym_imag[i];
 if(!isfinite(w)||w<=0)return 4;
 double measure=sqrt(3.129*(9.3-3.129))*mev;
 if(fabs(a.phase_space_J-measure)>1e-12*measure)return 5;
 puts("installed spectrum C11 caller passed");return 0;
}
