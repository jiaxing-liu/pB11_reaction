#include "fusion_boundary_segment.h"
#include <iostream>
#include <stdexcept>
void ck(bool b){if(!b)throw std::runtime_error("segment check failed");}
int main(){try{
 double coeff[]={2,.4,.4},a[]={2,0,0},b[]={2.1,0,.1};fusion_boundary_segment_value_v1 out{};
 auto run=[&](int N,double tol,int budget){return fusion_c_boundary_segment(1,coeff,N,a,b,1e-12,1e-12,tol,budget,&out);};
 ck(run(256,.001,2000)==0&&out.result==FUSION_SEGMENT_CLEAR);
 b[0]=2.8;b[2]=0;ck(run(256,1e-5,2000)==0&&out.result==FUSION_SEGMENT_GEOMETRIC_CANDIDATE);
 ck(out.fraction_lo<.5&&out.fraction_hi>.5); // exact torus exit at x=2.4
 ck(out.tolerance_met==0); // midpoint at the actual boundary is ambiguous
 ck(run(256,.001,2000)==0&&out.result==FUSION_SEGMENT_GEOMETRIC_CANDIDATE&&out.tolerance_met==0);
 ck(run(256,.0015,2000)==0&&out.result==FUSION_SEGMENT_GEOMETRIC_CANDIDATE&&out.tolerance_met==1);
 ck(out.fraction_lo<.5&&out.fraction_hi>.5&&out.fraction_hi-out.fraction_lo<=.0015);
 ck(run(256,1e-5,1)==0&&out.result==FUSION_SEGMENT_GEOMETRIC_CANDIDATE&&out.tolerance_met==0&&out.intervals_visited==1);
 // Both endpoints inside torus; segment crosses central hole twice.
 b[0]=-2;ck(run(256,.001,2000)==0&&out.result==FUSION_SEGMENT_GEOMETRIC_CANDIDATE);
 ck(out.fraction_lo<.1&&out.fraction_hi>.1&&out.fraction_hi<.2);
 // Endpoint on boundary, no safe outside witness: unresolved.
 b[0]=2.4;ck(run(256,1e-5,2000)==0&&out.result==FUSION_SEGMENT_UNRESOLVED);
 b[0]=-2;ck(run(256,1e-5,1)==0&&out.result==FUSION_SEGMENT_UNRESOLVED&&out.termination_reason==FUSION_SEGMENT_BUDGET_LIMIT);
 a[0]=2.4;ck(run(256,.001,2000)==0&&out.result==FUSION_SEGMENT_UNRESOLVED&&out.termination_reason==FUSION_SEGMENT_INITIAL_AMBIGUOUS);
 a[0]=3;ck(run(256,.001,2000)==PB11_STATUS_OUT_OF_RANGE&&out.fraction_hi==0);
 a[0]=2;ck(run(256,0,2000)==PB11_STATUS_INVALID_ARGUMENT);
 std::cout<<"segment inside, crossing, out-and-back, ambiguity, budget and errors passed\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
