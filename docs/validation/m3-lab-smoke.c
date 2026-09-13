#include "fusion_laboratory.h"
#include <math.h>
#include <stdio.h>
int main(void) {
 const double d[3]={0,0,1},v[3]={1e6,0,0};
 fusion_particle_four_vector_v1 p[2],q[2]; fusion_boost_ledger_v1 l;
 const double A=2e-12,ma=6.64e-27,mb=1.67e-27,c=299792458.;
 if(fusion_c_two_body_cm(ma,mb,A,d,p)) return 1;
 if(fabs(p[0].kinetic_energy_J+p[1].kinetic_energy_J-A)>1e-12*A) return 2;
 if(fusion_c_boost_particles(2,v,p,q,&l)) return 3;
 double g=1/sqrt(1-v[0]*v[0]/(c*c));
 double expected=A+(g*g*v[0]*v[0]/(c*c)/(g+1))*((ma+mb)*c*c+A);
 if(!isfinite(l.output_kinetic_energy_J) || fabs(l.output_kinetic_energy_J-expected)>1e-12*expected) return 4;
 if(fabs(l.energy_residual_J)>1e-12*expected) return 5;
 puts("installed laboratory C11 caller passed");return 0;
}
