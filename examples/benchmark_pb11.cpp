#include "pb11.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

namespace {

double relative_error(double fast, double integral) {
    return (fast - integral) / integral;
}

void print_energy_scan() {
    std::cout << "E_MeV,S_MeV_barn,sigma_barn\n";
    std::cout << std::scientific << std::setprecision(10);
    for (int energy_keV = 1; energy_keV <= 9760; ++energy_keV) {
        const double energy_MeV = energy_keV / 1000.0;
        std::cout << energy_MeV << ',' << pb11_sfactor(energy_MeV) << ','
                  << pb11_cross_section(energy_MeV) << '\n';
    }
}

void update_maximum(double temperature, double& max_error,
                    double& max_temperature) {
    const double integral = pb11_reactivity_integral(temperature);
    const double fast = pb11_reactivity_fast(temperature);
    const double absolute_error = std::abs(relative_error(fast, integral));
    if (absolute_error > max_error) {
        max_error = absolute_error;
        max_temperature = temperature;
    }
}

void print_reactivity_comparison() {
    constexpr std::array<double, 14> key_temperatures{
        10.0, 20.0, 30.0, 50.0, 60.0, 69.0, 69.9,
        70.0, 70.1, 100.0, 200.0, 300.0, 400.0, 500.0};

    std::cout << "T_keV,integral_m3_per_s,fast_m3_per_s,relative_error,"
                 "absolute_relative_error\n";
    std::cout << std::scientific << std::setprecision(10);
    for (double temperature : key_temperatures) {
        const double integral = pb11_reactivity_integral(temperature);
        const double fast = pb11_reactivity_fast(temperature);
        const double error = relative_error(fast, integral);
        std::cout << temperature << ',' << integral << ',' << fast << ','
                  << error << ',' << std::abs(error) << '\n';
    }

    double low_max_error = -1.0;
    double low_max_temperature = 0.0;
    for (int temperature_keV = 10; temperature_keV < 70; ++temperature_keV) {
        update_maximum(static_cast<double>(temperature_keV), low_max_error,
                       low_max_temperature);
    }
    update_maximum(69.9, low_max_error, low_max_temperature);

    double high_max_error = -1.0;
    double high_max_temperature = 0.0;
    for (int temperature_keV = 70; temperature_keV <= 500;
         temperature_keV += 10) {
        update_maximum(static_cast<double>(temperature_keV), high_max_error,
                       high_max_temperature);
    }
    update_maximum(70.1, high_max_error, high_max_temperature);

    std::cout << "# low_T_max_abs_error=" << low_max_error
              << ",T_keV=" << low_max_temperature << '\n';
    std::cout << "# high_T_max_abs_error=" << high_max_error
              << ",T_keV=" << high_max_temperature << '\n';
}

void print_reactivity_scan() {
    std::cout << "T_keV,integral_m3_per_s,fast_m3_per_s,relative_error,"
                 "absolute_relative_error\n";
    std::cout << std::scientific << std::setprecision(10);

    // A 1 keV grid is sufficient for the full curve.  Add 69.9 and 70.1 keV
    // explicitly to resolve the LT/HT switch around 70 keV.
    for (int temperature_tenths_keV = 100;
         temperature_tenths_keV <= 5000; ++temperature_tenths_keV) {
        if (temperature_tenths_keV % 10 != 0 &&
            temperature_tenths_keV != 699 &&
            temperature_tenths_keV != 701) {
            continue;
        }
        const double temperature = temperature_tenths_keV / 10.0;
        const double integral = pb11_reactivity_integral(temperature);
        const double fast = pb11_reactivity_fast(temperature);
        const double error = relative_error(fast, integral);
        std::cout << temperature << ',' << integral << ',' << fast << ','
                  << error << ',' << std::abs(error) << '\n';
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--energy-scan") {
        print_energy_scan();
        return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "--reactivity-scan") {
        print_reactivity_scan();
        return 0;
    }
    if (argc != 1) {
        std::cerr << "usage: benchmark_pb11 "
                     "[--energy-scan|--reactivity-scan]\n";
        return 2;
    }
    print_reactivity_comparison();
    return 0;
}
