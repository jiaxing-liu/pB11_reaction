#include "fusion_cross_sections.h"
#include "fusion_rates.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
constexpr double kev=1.602176634e-16;
int failures=0;
void check(bool condition,const char *name) {
    if (!condition) {++failures; std::cerr<<"FAIL: "<<name<<'\n';}
}
bool close(double a,double b,double tol) {return std::abs(a-b)<=tol*std::max(std::abs(a),std::abs(b));}
double integrated_rate(int channel,double temperature) {
    const double lower=channel==4?.3:.5;
    const double upper[]={0,5000,4900,4700,4800};
    const double mass[]={0,937814,937814,1124656,1124572};
    const double mu=mass[channel]*kev/(299792458.*299792458.);
    std::vector<double> points{lower,1,2,5,10,20,50,100,200,530,900,upper[channel]};
    std::sort(points.begin(),points.end());
    double integral=0;
    auto f=[=](double e) {
        double sigma;
        if (fusion_c_cross_section(channel,e*kev,&sigma)) throw std::runtime_error("integration data domain");
        return e*(sigma*1e28)*std::exp(-e/temperature);
    };
    for (std::size_t i=1;i<points.size();++i)
        integral+=boost::math::quadrature::gauss_kronrod<double,61>::integrate(f,points[i-1],points[i],15,1e-10);
    return std::sqrt(8*kev/(std::acos(-1.)*mu))*integral/std::pow(temperature,1.5)*1e-28;
}
}
int main() {
    // Table V p621: millibarn values, columns DT,DHe3,DD(Tp),DD(He3n).
    constexpr double energy[]={10,50,100};
    constexpr int channels[]={3,4,1,2};
    constexpr double expected[][4]={{27.02,2.160e-4,.2812,.2779},{4219,8.688,15.57,16.49},{3427,102.1,33.04,37.01}};
    for (int row=0;row<3;++row) for (int col=0;col<4;++col) {
        double sigma=-1;
        check(fusion_c_cross_section(channels[col],energy[row]*kev,&sigma)==0,"cross-section reference status");
        check(close(sigma,expected[row][col]*1e-31,6e-4),"Table V CM cross section and mb conversion");
    }
    for (int ch=1;ch<5;++ch) {
        double max_error=0;
        for (double temperature : {1.,5.,10.,20.,50.,100.}) {
            double fitted;
            check(fusion_c_thermal_reactivity(ch,temperature*kev,0,&fitted)==0,"integral comparison fitted status");
            const double numerical=integrated_rate(ch,temperature);
            max_error=std::max(max_error,std::abs(numerical/fitted-1));
            // Separate fitted cross sections and fitted rates have independent
            // Table IV/VI/VII errors. This is not a quadrature tolerance.
            check(close(numerical,fitted,ch==4?.05:.03),"cross-section quadrature versus thermal fit");
        }
        std::cout<<"channel "<<ch<<" maximum integral/fit relative difference "<<max_error<<'\n';
        double sigma=-1;
        check(fusion_c_cross_section(ch,0,&sigma)==0 && sigma==0,"zero energy limit");
        check(fusion_c_cross_section(ch,.1*kev,&sigma)==PB11_STATUS_OUT_OF_RANGE && sigma==0,"below positive data domain");
        check(fusion_c_cross_section(ch,6000*kev,&sigma)==PB11_STATUS_OUT_OF_RANGE && sigma==0,"above data domain");
    }
    // The DT fit should match near its documented 530 keV transition.
    double left,right;
    check(fusion_c_cross_section(3,529.999*kev,&left)==0 && fusion_c_cross_section(3,530.001*kev,&right)==0,"DT transition status");
    check(close(left,right,.002),"DT high fit uses constant numerator and matches low fit");
    check(fusion_c_cross_section(4,899.999*kev,&left)==0 && fusion_c_cross_section(4,900.001*kev,&right)==0,"DHe3 transition status");
    check(right>left && right/left<1.05,"DHe3 documented positive fit jump");
    check(fusion_c_cross_section(-1,kev,&left)==PB11_STATUS_INVALID_ARGUMENT && left==0,"bad channel clears cross section");
    check(fusion_c_cross_section(1,kev,nullptr)==PB11_STATUS_NULL_OUTPUT,"null cross section pointer");
    if (failures) return 1;
    std::cout<<"PASS: published cross sections, thermal quadrature, fit transitions and domains\n";
}
