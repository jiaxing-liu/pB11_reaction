#include "fusion_products.h"

#include <math.h>
#include <stdio.h>

static int close_to(double a, double b, double scale)
{
    return isfinite(a) && isfinite(b) && fabs(a - b) <= 1.0e-11 * fmax(1.0e-300, fabs(scale));
}

static int close_zero(double value, double scale)
{
    return isfinite(value) && fabs(value) <= 1.0e-9 * fmax(1.0e-30, fabs(scale));
}

int main(void)
{
    const double edges[3] = {0.0, 2.0e-15, 4.0e-15};
    const double packet_energy[2] = {2.0e-15, 5.0e-15};
    const double packet_rate[2] = {2.0e20, 3.0e20};
    double cell_birth[2] = {-1.0, -1.0};
    fusion_birth_mapping_v1 mapping;
    int status = fusion_c_map_birth_packets(2, edges, 2, packet_energy,
                                            packet_rate, cell_birth, &mapping);
    if (status != PB11_STATUS_OK ||
        !close_to(mapping.input_number_m3_s, 5.0e20, 5.0e20) ||
        !close_to(mapping.input_energy_W_m3, 1.9e6, 1.9e6) ||
        !close_to(mapping.mapped_number_m3_s, 2.0e20, 2.0e20) ||
        !close_to(mapping.mapped_energy_W_m3, 4.0e5, 4.0e5) ||
        !close_to(mapping.above_number_m3_s, 3.0e20, 3.0e20) ||
        !close_to(mapping.above_energy_W_m3, 1.5e6, 1.5e6) ||
        !close_to(cell_birth[0], 1.0e20, 1.0e20) ||
        !close_to(cell_birth[1], 1.0e20, 1.0e20) ||
        !close_zero(mapping.below_number_m3_s, 1.0) ||
        !close_zero(mapping.below_energy_W_m3, 1.0) ||
        !close_zero(mapping.number_residual_m3_s, 5.0e20) ||
        !close_zero(mapping.energy_residual_W_m3, 1.9e6)) {
        fprintf(stderr, "map smoke failed: status=%d mapped=%g above=%g birth=(%g,%g)\n",
                status, mapping.mapped_number_m3_s, mapping.above_number_m3_s,
                cell_birth[0], cell_birth[1]);
        return 1;
    }

    fusion_three_body_cm_v1 three_body;
    status = fusion_c_three_equal_sequential_cm(3.0e-27, 1.0e-13,
                                                 2.0e-14, -0.25, &three_body);
    double energy_sum = three_body.kinetic_energy_J[0]
                      + three_body.kinetic_energy_J[1]
                      + three_body.kinetic_energy_J[2];
    double px_sum = three_body.momentum_x_kg_m_s[0]
                  + three_body.momentum_x_kg_m_s[1]
                  + three_body.momentum_x_kg_m_s[2];
    double pz_sum = three_body.momentum_z_kg_m_s[0]
                  + three_body.momentum_z_kg_m_s[1]
                  + three_body.momentum_z_kg_m_s[2];
    double pscale = fabs(three_body.momentum_x_kg_m_s[0])
                  + fabs(three_body.momentum_x_kg_m_s[1])
                  + fabs(three_body.momentum_x_kg_m_s[2])
                  + fabs(three_body.momentum_z_kg_m_s[0])
                  + fabs(three_body.momentum_z_kg_m_s[1])
                  + fabs(three_body.momentum_z_kg_m_s[2]);
    if (status != PB11_STATUS_OK ||
        !close_to(energy_sum, 1.0e-13, 1.0e-13) ||
        !close_zero(px_sum, pscale) || !close_zero(pz_sum, pscale) ||
        !close_zero(three_body.energy_residual_J, 1.0e-13) ||
        !close_zero(three_body.momentum_residual_kg_m_s, pscale) ||
        !isfinite(three_body.kinetic_energy_J[0]) ||
        !isfinite(three_body.kinetic_energy_J[1]) ||
        !isfinite(three_body.kinetic_energy_J[2])) {
        fprintf(stderr, "three-body smoke failed: status=%d E=%g px=%g pz=%g\n",
                status, energy_sum, px_sum, pz_sum);
        return 2;
    }

    puts("fusion_products installed C11 smoke passed");
    return 0;
}
