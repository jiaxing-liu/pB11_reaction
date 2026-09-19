#include "fusion_boundary_segment.h"
#include <algorithm>
#include <cmath>
#include <limits>
namespace {
struct Point {double fraction,x[3];fusion_flux_boundary_value_v1 boundary;};
struct Node {Point a,b;int depth;};
}
namespace fusion_detail { double boundary_geometry_allowance(const fusion_flux_boundary_prepared*); }
template<class Query> int segment_impl(Query query,double geom,
 const double*start,const double*end,double tolerance,
 int budget,fusion_boundary_segment_value_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!start||!end||!std::isfinite(tolerance)||tolerance<=0||tolerance>1||budget<1||budget>1000000)return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<3;++j)if(!std::isfinite(start[j])||!std::isfinite(end[j]))return PB11_STATUS_INVALID_ARGUMENT;
 auto point=[&](double f,Point&p){
  p={};p.fraction=f;
  for(int j=0;j<3;++j){long double val=(1.L-f)*start[j]+static_cast<long double>(f)*end[j];
   if(!std::isfinite(val)||std::abs(val)>std::numeric_limits<double>::max())return int(PB11_STATUS_NUMERICAL_FAILURE);
   p.x[j]=double(val);
  }
  double R=std::hypot(p.x[0],p.x[1]);
  if(!std::isfinite(R))return int(PB11_STATUS_NUMERICAL_FAILURE);
  return query(R,p.x[2],&p.boundary);
 };
 Point first{},last{};int rc=point(0,first);if(rc)return rc;rc=point(1,last);if(rc)return rc;
 if(first.boundary.classification==FUSION_BOUNDARY_OUTSIDE)return PB11_STATUS_OUT_OF_RANGE;
 if(first.boundary.classification==FUSION_BOUNDARY_AMBIGUOUS){*out={0,1,first.boundary.uncertainty_m,FUSION_SEGMENT_UNRESOLVED,0,FUSION_SEGMENT_INITIAL_AMBIGUOUS,0};return 0;}
 Node stack[54];int size=1,visited=0;stack[0]={first,last,0};
 bool have_candidate=false;Node candidate{};
 auto finish=[&](const Node&leaf,int reason){
  if(have_candidate){
   // leaf.a is reached only after every preceding interval was covered inside.
   if(leaf.a.boundary.classification!=FUSION_BOUNDARY_INSIDE||leaf.a.fraction<candidate.a.fraction||leaf.a.fraction>=candidate.b.fraction)return int(PB11_STATUS_NUMERICAL_FAILURE);
   candidate.a=leaf.a;
   double search_lo=candidate.a.fraction;
   while(candidate.b.fraction-candidate.a.fraction>tolerance&&visited<budget){
    double f=search_lo+(candidate.b.fraction-search_lo)/2;
    if(!(f>search_lo&&f<candidate.b.fraction)){reason=FUSION_SEGMENT_RESOLUTION_LIMIT;break;}
    Point p{};int err=point(f,p);++visited;if(err)return err;
    if(p.boundary.classification==FUSION_BOUNDARY_OUTSIDE)candidate.b=p;
    else search_lo=f; // search hint only: NEVER discard the unproved prefix
   }
   if(candidate.b.fraction-candidate.a.fraction<=tolerance)reason=FUSION_SEGMENT_FRACTION_LIMIT;
   else if(visited>=budget)reason=FUSION_SEGMENT_BUDGET_LIMIT;
  }
  const Node&r=have_candidate?candidate:leaf;
  *out={r.a.fraction,r.b.fraction,first.boundary.uncertainty_m,
   have_candidate?FUSION_SEGMENT_GEOMETRIC_CANDIDATE:FUSION_SEGMENT_UNRESOLVED,
   (r.b.fraction-r.a.fraction)<=tolerance?1:0,reason,visited};
  return int(PB11_STATUS_OK);
 };
 while(size){
  if(visited>=budget)return finish(stack[size-1],FUSION_SEGMENT_BUDGET_LIMIT);
  Node node=stack[--size];++visited;
  const auto&a=node.a;const auto&b=node.b;
  bool inside_a=a.boundary.classification==FUSION_BOUNDARY_INSIDE;
  bool inside_b=b.boundary.classification==FUSION_BOUNDARY_INSIDE;
  if(!inside_a)return PB11_STATUS_NUMERICAL_FAILURE; // prior prefix was covered
  long double dx=static_cast<long double>(b.x[0])-a.x[0],dy=static_cast<long double>(b.x[1])-a.x[1],dz=static_cast<long double>(b.x[2])-a.x[2];
  long double len=std::hypot(std::hypot(dx,dy),dz)+geom;
  long double clearance=std::max(static_cast<long double>(a.boundary.polygon_distance_m)-a.boundary.uncertainty_m,static_cast<long double>(b.boundary.polygon_distance_m)-b.boundary.uncertainty_m);
  if(inside_b&&clearance>len)continue;
  if(b.boundary.classification==FUSION_BOUNDARY_OUTSIDE){candidate=node;have_candidate=true;}
  if(b.fraction-a.fraction<=tolerance)return finish(node,FUSION_SEGMENT_FRACTION_LIMIT);
  if(visited>=budget)return finish(node,FUSION_SEGMENT_BUDGET_LIMIT);
  if(node.depth>=52)return finish(node,FUSION_SEGMENT_RESOLUTION_LIMIT);
  double f=a.fraction+(b.fraction-a.fraction)/2;
  if(!(f>a.fraction&&f<b.fraction))return finish(node,FUSION_SEGMENT_RESOLUTION_LIMIT);
  Point middle{};rc=point(f,middle);if(rc)return rc;
  stack[size++]={middle,b,node.depth+1};stack[size++]={a,middle,node.depth+1};
 }
 *out={1,1,first.boundary.uncertainty_m,FUSION_SEGMENT_CLEAR,1,FUSION_SEGMENT_FINISHED,visited};return 0;
}

extern "C" int fusion_c_boundary_segment(int m,const double*coeff,int segments,
 const double*start,const double*end,double geom,double angular,double tolerance,
 int budget,fusion_boundary_segment_value_v1*out){
 auto query=[&](double R,double Z,fusion_flux_boundary_value_v1*b){return fusion_c_flux_boundary(m,coeff,segments,R,Z,geom,angular,b);};
 return segment_impl(query,geom,start,end,tolerance,budget,out);
}
extern "C" int fusion_c_boundary_prepared_segment(const fusion_flux_boundary_prepared*context,
 const double*start,const double*end,double tolerance,int budget,fusion_boundary_segment_value_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};if(!context)return PB11_STATUS_INVALID_ARGUMENT;
 auto query=[&](double R,double Z,fusion_flux_boundary_value_v1*b){return fusion_c_flux_boundary_prepared_point(context,R,Z,b);};
 return segment_impl(query,fusion_detail::boundary_geometry_allowance(context),start,end,tolerance,budget,out);
}
