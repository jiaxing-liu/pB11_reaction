#include "fusion_flux_boundary.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void ck(bool b){if(!b)throw std::runtime_error("boundary check failed");}
int main(){try{
 const double pi=std::acos(-1.);fusion_flux_boundary_value_v1 o{};
 for(double b:{.4,.6}){
  double a[]={2,.4,b};double old=0;
  for(int N:{64,128,256}){
   for(int j=0;j<128;++j){double t=j*2*pi/128;
    for(double scale:{.95,1.,1.05}){
     ck(fusion_c_flux_boundary(1,a,N,2+.4*scale*std::cos(t),b*scale*std::sin(t),1e-12,1e-12,&o)==0);
     ck(o.classification==(scale<1?FUSION_BOUNDARY_INSIDE:scale>1?FUSION_BOUNDARY_OUTSIDE:FUSION_BOUNDARY_AMBIGUOUS));
    }
   }
   if(old)ck(std::abs((old-1e-12)/(o.uncertainty_m-1e-12)-4)<1e-12);
   old=o.uncertainty_m;
  }
 }
 // Shaped fifth harmonic, still within the supported angular-monotone domain.
 double high[]={2,.4,0,0,0,.003,.6,0,0,0,-.002};
 for(int j=0;j<1000;++j){double t=j*2*pi/1000;
  ck(fusion_c_flux_boundary(5,high,256,2+.4*std::cos(t)+.003*std::cos(5*t),.6*std::sin(t)-.002*std::sin(5*t),1e-12,1e-12,&o)==0);
  ck(o.classification==FUSION_BOUNDARY_AMBIGUOUS);
 }
 double twice[]={2,0,.4,0,.6};
 ck(fusion_c_flux_boundary(2,twice,256,2,0,1e-12,1e-12,&o)==FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY);
 ck(o.classification==0&&o.uncertainty_m==0);
 double reverse[]={2,.4,-.6};
 ck(fusion_c_flux_boundary(1,reverse,256,2,0,1e-12,1e-12,&o)==FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY);
 double degenerate[]={2,.4,0};
 ck(fusion_c_flux_boundary(1,degenerate,256,2,0,1e-12,1e-12,&o)==FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY);
 ck(fusion_c_flux_boundary(1,reverse,256,2,0,0,1e-12,&o)==PB11_STATUS_INVALID_ARGUMENT);
 ck(fusion_c_flux_boundary(1,nullptr,256,2,0,1e-12,1e-12,&o)==PB11_STATUS_INVALID_ARGUMENT);
 ck(fusion_c_flux_boundary(1,reverse,256,2,0,1e-12,1e-12,nullptr)==PB11_STATUS_NULL_OUTPUT);
 std::cout<<"boundary analytic, shaped, refinement, winding and failure checks passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
