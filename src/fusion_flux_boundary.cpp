#include "fusion_flux_boundary.h"
#include <algorithm>
#include <cmath>
#include <limits>
namespace {
using L=long double;
struct Point {L r,z,dr,dz;};
Point point(int m,const double*a,L t){
 Point p{};
 for(int j=1;j<=m;++j){L co=std::cos(j*t),si=std::sin(j*t);
  p.r+=a[j]*co;p.z+=a[m+j]*si;
  p.dr-=j*L(a[j])*si;p.dz+=j*L(a[m+j])*co;
 }
 return p;
}
L distance(Point a,Point b,L R,L Z){
 L u=b.r-a.r,v=b.z-a.z,den=u*u+v*v;
 if(!(den>0))return -1;
 L t=std::max(0.L,std::min(1.L,((R-a.r)*u+(Z-a.z)*v)/den));
 return std::hypot(R-a.r-t*u,Z-a.z-t*v);
}
int crossing(Point a,Point b,L R,L Z){
 L cross=(b.r-a.r)*(Z-a.z)-(b.z-a.z)*(R-a.r);
 if(a.z<=Z&&b.z>Z&&cross>0)return 1;
 if(a.z>Z&&b.z<=Z&&cross<0)return -1;
 return 0;
}
}
extern "C" int fusion_c_flux_boundary(int m,const double*a,int segments,
 double R,double Z,double geom,double angular,fusion_flux_boundary_value_v1*out){
 if(!out)return PB11_STATUS_NULL_OUTPUT;
 *out={};
 if(!a||m<1||m>64||segments<16||segments>65536||!std::isfinite(R)||!std::isfinite(Z)||!std::isfinite(geom)||geom<=0||!std::isfinite(angular)||angular<=0)return PB11_STATUS_INVALID_ARGUMENT;
 for(int j=0;j<=2*m;++j)if(!std::isfinite(a[j]))return PB11_STATUS_INVALID_ARGUMENT;
 L P=0,M=0,rbound=0;
 for(int j=1;j<=m;++j){L amp=std::max(std::abs(L(a[j])),std::abs(L(a[m+j])));P+=amp;M+=j*j*amp;rbound+=std::abs(L(a[j]));}
 L h=2*std::acos(-1.L)/segments,band=M*h*h/8+geom;
 if(!(L(a[0])-rbound>geom))return FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY;
 L minW=std::numeric_limits<L>::infinity(),dc=minW,dq=minW;
 int wc=0,wq=0;L qr=L(R)-a[0],qz=Z;
 Point first=point(m,a,0),prev=first;
 for(int i=0;i<segments;++i){
  Point next=i+1==segments?first:point(m,a,(i+1)*h);
  minW=std::min(minW,prev.r*prev.dz-prev.z*prev.dr);
  L cdist=distance(prev,next,0,0),qdist=distance(prev,next,qr,qz);
  if(cdist<0||qdist<0)return FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY;
  dc=std::min(dc,cdist);dq=std::min(dq,qdist);
  wc+=crossing(prev,next,0,0);wq+=crossing(prev,next,qr,qz);prev=next;
 }
 L margin=minW-P*M*h/2-angular;
 for(L val:{band,margin,dc,dq})if(!std::isfinite(val)||std::abs(val)>std::numeric_limits<double>::max())return PB11_STATUS_NUMERICAL_FAILURE;
 if(!(margin>0)||!(dc>band)||wc!=1)return FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY;
 int classification=FUSION_BOUNDARY_AMBIGUOUS;
 if(dq>band){
  if(wq!=0&&wq!=1)return FUSION_BOUNDARY_UNSUPPORTED_GEOMETRY;
  classification=wq==1?FUSION_BOUNDARY_INSIDE:FUSION_BOUNDARY_OUTSIDE;
 }
 *out={double(dq),double(band),double(margin),classification};
 return PB11_STATUS_OK;
}
