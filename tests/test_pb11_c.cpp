#include "pb11_c.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

}  // namespace

int main() {
    check(pb11_c_abi_version() == 1, "C ABI version is 1");

    double value = -1.0;
    int status = pb11_c_sfactor(0.148, &value);
    check(status == PB11_STATUS_OK, "valid S-factor call returns OK");
    check(std::isfinite(value) && value > 3000.0,
          "valid S-factor call returns a finite value");

    status = pb11_c_cross_section(0.0, &value);
    check(status == PB11_STATUS_OK && value == 0.0,
          "zero-energy cross section is accepted through the C ABI");

    value = -1.0;
    status = pb11_c_reactivity_integral(100.0, &value);
    check(status == PB11_STATUS_OK && std::isfinite(value) && value > 0.0,
          "integral reactivity call returns a finite positive value");

    double generic_integral = 0.0;
    status = pb11_c_reactivity(100.0, PB11_REACTIVITY_INTEGRAL,
                               &generic_integral);
    check(status == PB11_STATUS_OK && generic_integral == value,
          "generic integral method selects the integral implementation");

    status = pb11_c_reactivity_fast(100.0, &value);
    check(status == PB11_STATUS_OK && std::isfinite(value) && value > 0.0,
          "fast reactivity call returns a finite positive value");

    value = -1.0;
    status = pb11_c_reactivity_fast(9.999, &value);
    check(status == PB11_STATUS_OUT_OF_RANGE && value == 0.0,
          "fast method rejects temperatures below its published range");

    value = -1.0;
    status = pb11_c_reactivity_integral(0.0, &value);
    check(status == PB11_STATUS_OUT_OF_RANGE && value == 0.0,
          "integral method rejects non-positive temperatures");

    value = -1.0;
    status = pb11_c_sfactor(std::numeric_limits<double>::quiet_NaN(), &value);
    check(status == PB11_STATUS_INVALID_ARGUMENT && value == 0.0,
          "non-finite input is rejected without returning NaN");

    status = pb11_c_sfactor(0.1, nullptr);
    check(status == PB11_STATUS_NULL_OUTPUT,
          "null output pointer is reported");

    value = -1.0;
    status = pb11_c_reactivity(100.0, 99, &value);
    check(status == PB11_STATUS_UNKNOWN_METHOD && value == 0.0,
          "unknown reactivity method is reported");

    check(std::string(pb11_c_status_message(PB11_STATUS_OUT_OF_RANGE)) ==
              "input is outside the published range",
          "status message is stable for range errors");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All p-11B C ABI tests passed\n";
    return 0;
}
