program test_fusion_rate_model_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_beam_fortran, only : fusion_beam_window_v1, FUSION_DT_ALPHAN
  use fusion_rate_model_fortran
  implicit none

  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  real(c_double), parameter :: deuteron_mass = 3.3435837724e-27_c_double
  real(c_double), parameter :: triton_mass = 5.0073567446e-27_c_double
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_cross_section_model()
  call verify_thermal_model_and_segment_sum()
  call verify_cold_beam_below_domain()
  call verify_invalid_policy_clears_outputs()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran fusion-rate-model test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran fusion-rate-model tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
       write(*, '(A)') 'FAIL: ' // message
       failures = failures + 1
    end if
  end subroutine check

  logical function close_relative(actual, expected, tolerance)
    real(c_double), intent(in) :: actual, expected, tolerance
    real(c_double) :: scale

    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
       close_relative = .false.
       return
    end if
    if (actual == expected) then
       close_relative = .true.
       return
    end if
    scale = max(abs(actual), abs(expected))
    close_relative = scale > 0.0_c_double .and. &
         abs(actual - expected) <= tolerance * scale
  end function close_relative

  logical function window_is_zero(out)
    type(fusion_beam_window_v1), intent(in) :: out

    window_is_zero = out%resolved_reactivity_m3_s == 0.0_c_double .and. &
         out%projectile_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%target_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%relative_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%cm_energy_reactivity_J_m3_s == 0.0_c_double .and. &
         out%resolved_pair_probability == 0.0_c_double .and. &
         out%unresolved_pair_probability == 0.0_c_double .and. &
         out%unresolved_relative_speed_m_s == 0.0_c_double .and. &
         out%quadrature_error_m3_s == 0.0_c_double .and. &
         out%energy_identity_error_J_m3_s == 0.0_c_double .and. &
         out%domain_incomplete == 0_c_int
  end function window_is_zero

  logical function model_is_zero(out)
    type(fusion_rate_model_v1), intent(in) :: out

    model_is_zero = window_is_zero(out%total) .and. &
         window_is_zero(out%fit) .and. window_is_zero(out%below) .and. &
         window_is_zero(out%above)
  end function model_is_zero

  logical function window_is_finite(out)
    type(fusion_beam_window_v1), intent(in) :: out

    window_is_finite = ieee_is_finite(out%resolved_reactivity_m3_s) .and. &
         ieee_is_finite(out%projectile_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%target_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%relative_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%cm_energy_reactivity_J_m3_s) .and. &
         ieee_is_finite(out%resolved_pair_probability) .and. &
         ieee_is_finite(out%unresolved_pair_probability) .and. &
         ieee_is_finite(out%unresolved_relative_speed_m_s) .and. &
         ieee_is_finite(out%quadrature_error_m3_s) .and. &
         ieee_is_finite(out%energy_identity_error_J_m3_s)
  end function window_is_finite

  subroutine verify_layout()
    type(fusion_beam_window_v1) :: window
    type(fusion_rate_model_v1) :: model

    call check(c_sizeof(model) == 4_c_size_t * c_sizeof(window), &
         'rate-model struct has four beam-window values')
  end subroutine verify_layout

  subroutine verify_cross_section_model()
    real(c_double) :: energy, sigma
    integer(c_int) :: status

    energy = 10.0_c_double * joules_per_keV
    sigma = -1.0_c_double
    call fusion_cross_section_model(FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, &
         FUSION_PB_LOW_TB, energy, sigma, status)
    call check(status == PB11_STATUS_OK .and. sigma > 0.0_c_double .and. &
         ieee_is_finite(sigma), 'cross-section model wrapper returns a value')

    sigma = -1.0_c_double
    call fusion_cross_section_model(FUSION_DT_ALPHAN, 0_c_int, &
         FUSION_PB_LOW_TB, energy, sigma, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. &
         sigma == 0.0_c_double, &
         'invalid cross-section policy clears its scalar output')
  end subroutine verify_cross_section_model

  subroutine verify_thermal_model_and_segment_sum()
    real(c_double), parameter :: temperature = 20.0_c_double * &
         joules_per_keV
    type(fusion_rate_model_v1) :: out
    integer(c_int) :: status

    call fusion_thermal_pair_maxwellian_model(FUSION_DT_ALPHAN, &
         FUSION_ENDPOINT_S, FUSION_PB_LOW_TB, deuteron_mass, triton_mass, &
         temperature, temperature, out, status)
    call check(status == PB11_STATUS_OK, &
         'thermal DT rate-model wrapper returns OK')
    call check(window_is_finite(out%total) .and. window_is_finite(out%fit) .and. &
         window_is_finite(out%below) .and. window_is_finite(out%above), &
         'thermal DT model windows are finite')
    call check(out%total%resolved_reactivity_m3_s > 0.0_c_double, &
         'thermal DT total reactivity is positive')
    call check(out%below%resolved_reactivity_m3_s > 0.0_c_double, &
         'thermal DT has a nonzero below-domain tail')
    call check(out%total%domain_incomplete == 0_c_int, &
         'bounded model total has a complete continuation domain')
    call check_model_segment_sums(out)
  end subroutine verify_thermal_model_and_segment_sum

  subroutine check_model_segment_sums(out)
    type(fusion_rate_model_v1), intent(in) :: out
    real(c_double) :: expected

    expected = out%fit%resolved_reactivity_m3_s + &
         out%below%resolved_reactivity_m3_s + &
         out%above%resolved_reactivity_m3_s
    call check(close_relative(out%total%resolved_reactivity_m3_s, expected, &
         1.0e-12_c_double), 'total reactivity equals segment sum')

    expected = out%fit%projectile_energy_reactivity_J_m3_s + &
         out%below%projectile_energy_reactivity_J_m3_s + &
         out%above%projectile_energy_reactivity_J_m3_s
    call check(close_relative(out%total%projectile_energy_reactivity_J_m3_s, &
         expected, 1.0e-12_c_double), 'total projectile moment equals sum')

    expected = out%fit%target_energy_reactivity_J_m3_s + &
         out%below%target_energy_reactivity_J_m3_s + &
         out%above%target_energy_reactivity_J_m3_s
    call check(close_relative(out%total%target_energy_reactivity_J_m3_s, &
         expected, 1.0e-12_c_double), 'total target moment equals sum')

    expected = out%fit%relative_energy_reactivity_J_m3_s + &
         out%below%relative_energy_reactivity_J_m3_s + &
         out%above%relative_energy_reactivity_J_m3_s
    call check(close_relative(out%total%relative_energy_reactivity_J_m3_s, &
         expected, 1.0e-12_c_double), 'total relative moment equals sum')

    expected = out%fit%cm_energy_reactivity_J_m3_s + &
         out%below%cm_energy_reactivity_J_m3_s + &
         out%above%cm_energy_reactivity_J_m3_s
    call check(close_relative(out%total%cm_energy_reactivity_J_m3_s, expected, &
         1.0e-12_c_double), 'total CM moment equals sum')

    expected = out%fit%resolved_pair_probability + &
         out%below%resolved_pair_probability + &
         out%above%resolved_pair_probability
    call check(close_relative(out%total%resolved_pair_probability, expected, &
         1.0e-12_c_double), 'total resolved probability equals segment sum')

    expected = out%fit%quadrature_error_m3_s + &
         out%below%quadrature_error_m3_s + &
         out%above%quadrature_error_m3_s
    call check(close_relative(out%total%quadrature_error_m3_s, expected, &
         1.0e-12_c_double), 'total quadrature error equals segment sum')
  end subroutine check_model_segment_sums

  subroutine verify_cold_beam_below_domain()
    real(c_double), parameter :: relative_energy = 0.25_c_double * &
         joules_per_keV
    real(c_double) :: projectile_energy
    type(fusion_rate_model_v1) :: out
    integer(c_int) :: status

    projectile_energy = relative_energy * &
         (deuteron_mass + triton_mass) / triton_mass
    call fusion_beam_maxwellian_model(FUSION_DT_ALPHAN, FUSION_ENDPOINT_S, &
         FUSION_PB_LOW_TB, deuteron_mass, triton_mass, projectile_energy, &
         0.0_c_double, out, status)
    call check(status == PB11_STATUS_OK, &
         'cold beam below-domain model call returns OK')
    call check(out%below%resolved_reactivity_m3_s > 0.0_c_double .and. &
         out%fit%resolved_reactivity_m3_s == 0.0_c_double .and. &
         out%above%resolved_reactivity_m3_s == 0.0_c_double, &
         'cold beam below-domain rate is assigned to the below segment')
    call check(out%total%resolved_reactivity_m3_s == &
         out%below%resolved_reactivity_m3_s, &
         'cold beam total is the below-domain segment')
    call check(out%total%resolved_pair_probability == 1.0_c_double, &
         'cold beam model resolves the below-domain pair')
  end subroutine verify_cold_beam_below_domain

  subroutine verify_invalid_policy_clears_outputs()
    real(c_double), parameter :: temperature = 20.0_c_double * &
         joules_per_keV
    type(fusion_rate_model_v1) :: out
    integer(c_int) :: status

    call fusion_thermal_pair_maxwellian_model(FUSION_DT_ALPHAN, &
         FUSION_ENDPOINT_S, FUSION_PB_LOW_TB, deuteron_mass, triton_mass, &
         temperature, temperature, out, status)
    call check(status == PB11_STATUS_OK .and. &
         out%total%resolved_reactivity_m3_s > 0.0_c_double, &
         'valid model call populates output before error check')

    call fusion_thermal_pair_maxwellian_model(FUSION_DT_ALPHAN, 0_c_int, &
         FUSION_PB_LOW_TB, deuteron_mass, triton_mass, temperature, &
         temperature, out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT .and. model_is_zero(out), &
         'invalid model policy clears every output window')
  end subroutine verify_invalid_policy_clears_outputs

end program test_fusion_rate_model_fortran
