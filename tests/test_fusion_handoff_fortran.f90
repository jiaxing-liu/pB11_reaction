program test_fusion_handoff_fortran
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t, c_sizeof
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fusion_handoff_fortran
  implicit none

  real(c_double), parameter :: joules_per_keV = 1.602176634e-16_c_double
  real(c_double), parameter :: pi = 3.1415926535897932384626433832795_c_double
  integer :: failures

  failures = 0
  call verify_layout()
  call verify_grid_cdf()
  call verify_empty_handoff()
  call verify_one_cell_projection()
  call verify_bad_extents()

  if (failures /= 0) then
     write(*, '(I0, A)') failures, ' Fortran handoff test(s) failed'
     error stop 1
  end if
  write(*, '(A)') 'All Fortran handoff tests passed'

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

  logical function grid_is_zero(out)
    type(fusion_maxwellian_grid_v1), intent(in) :: out

    grid_is_zero = out%below_probability == 0.0_c_double .and. &
         out%above_probability == 0.0_c_double .and. &
         out%represented_mean_energy_J == 0.0_c_double .and. &
         out%probability_balance_error == 0.0_c_double
  end function grid_is_zero

  logical function handoff_is_zero(out)
    type(fusion_handoff_ledger_v1), intent(in) :: out

    handoff_is_zero = out%initial_number_m3 == 0.0_c_double .and. &
         out%initial_energy_J_m3 == 0.0_c_double .and. &
         out%remaining_number_m3 == 0.0_c_double .and. &
         out%remaining_energy_J_m3 == 0.0_c_double .and. &
         out%fluid_number_m3 == 0.0_c_double .and. &
         out%fluid_energy_J_m3 == 0.0_c_double .and. &
         out%bath_energy_correction_J_m3 == 0.0_c_double .and. &
         out%distribution_L1 == 0.0_c_double .and. &
         out%relative_mean_energy_error == 0.0_c_double .and. &
         out%outside_grid_probability == 0.0_c_double .and. &
         out%represented_Maxwellian_mean_energy_J == 0.0_c_double .and. &
         out%particle_balance_error_m3 == 0.0_c_double .and. &
         out%energy_balance_error_J_m3 == 0.0_c_double
  end function handoff_is_zero

  real(c_double) function maxwellian_cdf(x)
    real(c_double), intent(in) :: x

    if (x <= 0.0_c_double) then
       maxwellian_cdf = 0.0_c_double
    else
       maxwellian_cdf = erf(sqrt(x)) - 2.0_c_double * sqrt(x) * &
            exp(-x) / sqrt(pi)
    end if
  end function maxwellian_cdf

  real(c_double) function maxwellian_survival(x)
    real(c_double), intent(in) :: x

    if (x <= 0.0_c_double) then
       maxwellian_survival = 1.0_c_double
    else
       maxwellian_survival = erfc(sqrt(x)) + 2.0_c_double * sqrt(x) * &
            exp(-x) / sqrt(pi)
    end if
  end function maxwellian_survival

  subroutine verify_layout()
    type(fusion_maxwellian_grid_v1) :: grid
    type(fusion_handoff_ledger_v1) :: ledger
    real(c_double) :: d

    call check(c_sizeof(grid) == 4_c_size_t * c_sizeof(d), &
         'Maxwellian grid struct has four c_double fields')
    call check(c_sizeof(ledger) == 13_c_size_t * c_sizeof(d), &
         'handoff ledger struct has thirteen c_double fields')
  end subroutine verify_layout

  subroutine verify_grid_cdf()
    integer(c_int), parameter :: cells = 3_c_int
    real(c_double) :: edges(cells + 1), probability(cells)
    real(c_double) :: expected(cells), target, x0, x1, x2, x3
    real(c_double) :: total, represented
    type(fusion_maxwellian_grid_v1) :: out
    integer(c_int) :: status
    integer :: i

    target = 10.0_c_double * joules_per_keV
    edges = [0.0_c_double, 0.5_c_double * target, &
         1.5_c_double * target, 3.0_c_double * target]
    x0 = edges(1) / target
    x1 = edges(2) / target
    x2 = edges(3) / target
    x3 = edges(4) / target
    expected(1) = maxwellian_cdf(x1) - maxwellian_cdf(x0)
    expected(2) = maxwellian_cdf(x2) - maxwellian_cdf(x1)
    expected(3) = maxwellian_cdf(x3) - maxwellian_cdf(x2)
    probability = -1.0_c_double

    call fusion_maxwellian_energy_grid(cells, target, edges, probability, &
         out, status)
    call check(status == PB11_STATUS_OK, &
         'medium-x Maxwellian grid call returns OK')
    do i = 1, cells
       call check(close_relative(probability(i), expected(i), &
            2.0e-14_c_double), 'grid bin agrees with analytic 3D CDF')
    end do
    call check(close_relative(out%below_probability, &
         maxwellian_cdf(x0), 2.0e-14_c_double), 'grid below tail agrees with CDF')
    call check(close_relative(out%above_probability, &
         maxwellian_survival(x3), 2.0e-14_c_double), &
         'grid above tail agrees with survival CDF')
    total = sum(probability) + out%below_probability + out%above_probability
    represented = sum(probability * 0.5_c_double * &
         (edges(1:cells) + edges(2:cells + 1)))
    call check(close_relative(total, 1.0_c_double, 2.0e-14_c_double), &
         'grid CDF bins and tails close normalization')
    call check(close_relative(out%represented_mean_energy_J, represented, &
         2.0e-14_c_double), 'grid represented mean uses arithmetic centers')
    call check(abs(out%probability_balance_error) < 2.0e-14_c_double, &
         'grid probability residual is small')
    if (ieee_is_finite(out%represented_mean_energy_J)) then
       call check(.true., 'grid represented mean is finite')
    else
       call check(.false., 'grid represented mean is finite')
    end if
  end subroutine verify_grid_cdf

  subroutine verify_empty_handoff()
    real(c_double) :: edges(4), old_number(3), trial_number(3)
    type(fusion_handoff_ledger_v1) :: out
    integer(c_int) :: projected, status

    edges = [0.0_c_double, 0.5_c_double * joules_per_keV, &
         2.0_c_double * joules_per_keV, 5.0_c_double * joules_per_keV]
    old_number = 0.0_c_double
    trial_number = -1.0_c_double
    projected = -1_c_int

    call fusion_maxwellian_handoff_trial(3_c_int, 1.0_c_double * &
         joules_per_keV, 2.0_c_double, 0.01_c_double, edges, old_number, &
         trial_number, projected, out, status)
    call check(status == PB11_STATUS_OK, 'empty handoff trial returns OK')
    call check(projected == 0_c_int .and. all(trial_number == old_number), &
         'empty candidate is retained and never projected')
    call check(out%fluid_number_m3 == 0.0_c_double .and. &
         out%fluid_energy_J_m3 == 0.0_c_double .and. &
         out%bath_energy_correction_J_m3 == 0.0_c_double, &
         'empty handoff has no fluid source or bath correction')
    call check(out%remaining_number_m3 == 0.0_c_double .and. &
         out%remaining_energy_J_m3 == 0.0_c_double, &
         'empty handoff remaining source is zero')
  end subroutine verify_empty_handoff

  subroutine verify_one_cell_projection()
    real(c_double) :: target, edges(2), old_number(1), trial_number(1)
    real(c_double) :: old_before, initial_number, initial_energy
    real(c_double) :: second_trial(1)
    type(fusion_handoff_ledger_v1) :: out, second_out
    integer(c_int) :: projected, second_projected, status, second_status

    target = 10.0_c_double * joules_per_keV
    edges = [target, 2.0_c_double * target]
    old_number = [2.0e19_c_double]
    old_before = old_number(1)
    trial_number = -1.0_c_double
    projected = -1_c_int

    call fusion_maxwellian_handoff_trial(1_c_int, target, 2.0_c_double, &
         0.01_c_double, edges, old_number, trial_number, projected, out, &
         status)
    call check(status == PB11_STATUS_OK .and. projected == 1_c_int, &
         'one-cell 1.5-kT candidate is projected under explicit bounds')
    call check(trial_number(1) == 0.0_c_double, &
         'projected trial removes the kinetic candidate')

    initial_number = old_number(1)
    initial_energy = 1.5_c_double * target * initial_number
    call check(close_relative(out%initial_number_m3, initial_number, &
         1.0e-14_c_double), 'handoff ledger records source particle amount')
    call check(close_relative(out%initial_energy_J_m3, initial_energy, &
         1.0e-14_c_double), 'one-cell center gives 1.5-kT source energy')
    call check(close_relative(out%fluid_number_m3, initial_number, &
         1.0e-14_c_double), 'projected fluid particle source is Nold')
    call check(close_relative(out%fluid_energy_J_m3, initial_energy, &
         1.0e-14_c_double), 'projected fluid energy source is 1.5-kT Nold')
    call check(out%remaining_number_m3 == 0.0_c_double .and. &
         out%remaining_energy_J_m3 == 0.0_c_double, &
         'projected handoff has no remaining kinetic source')
    ! The C ledger balances actual rounded fluid energy against a long-double
    ! midpoint moment. Its correction need not be bitwise zero even here.
    call check(abs(out%bath_energy_correction_J_m3) <= &
         4.0_c_double * epsilon(1.0_c_double) * abs(initial_energy), &
         '1.5-kT cell correction is bounded by floating-point rounding')
    call check(close_relative(out%fluid_number_m3 + out%remaining_number_m3, &
         out%initial_number_m3, 1.0e-14_c_double), &
         'handoff particle source accounting closes')
    call check(close_relative(out%fluid_energy_J_m3 + &
         out%bath_energy_correction_J_m3, out%initial_energy_J_m3, &
         1.0e-14_c_double), 'handoff energy plus correction accounting closes')
    call check(abs(out%particle_balance_error_m3) < 1.0e-12_c_double .and. &
         abs(out%energy_balance_error_J_m3) < 1.0e-12_c_double, &
         'projected source ledger residuals are small')
    call check(old_number(1) == old_before, &
         'handoff projection leaves old candidate unchanged')

    call fusion_maxwellian_handoff_trial(1_c_int, target, 2.0_c_double, &
         0.01_c_double, edges, old_number, second_trial, second_projected, &
         second_out, second_status)
    call check(second_status == PB11_STATUS_OK .and. &
         second_projected == projected .and. second_trial(1) == trial_number(1), &
         'repeated handoff trial is deterministic')
    call check(second_out%initial_number_m3 == out%initial_number_m3 .and. &
         second_out%fluid_energy_J_m3 == out%fluid_energy_J_m3 .and. &
         old_number(1) == old_before, &
         'repeated handoff does not accumulate state or modify old input')
  end subroutine verify_one_cell_projection

  subroutine verify_bad_extents()
    real(c_double) :: grid_edges(3), short_probability(1)
    real(c_double) :: handoff_edges(3), old_number(2), short_trial(1)
    type(fusion_maxwellian_grid_v1) :: grid_out
    type(fusion_handoff_ledger_v1) :: handoff_out
    integer(c_int) :: projected, status

    grid_edges = [0.0_c_double, 1.0_c_double * joules_per_keV, &
         2.0_c_double * joules_per_keV]
    short_probability = -1.0_c_double
    call fusion_maxwellian_energy_grid(2_c_int, joules_per_keV, grid_edges, &
         short_probability, grid_out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'grid wrapper rejects a wrong probability extent')
    call check(short_probability(1) == 0.0_c_double .and. grid_is_zero(grid_out), &
         'grid extent failure clears every output')

    handoff_edges = grid_edges
    old_number = [1.0e19_c_double, 2.0e19_c_double]
    short_trial = -1.0_c_double
    projected = -1_c_int
    call fusion_maxwellian_handoff_trial(2_c_int, joules_per_keV, 2.0_c_double, &
         0.1_c_double, handoff_edges, old_number, short_trial, projected, &
         handoff_out, status)
    call check(status == PB11_STATUS_INVALID_ARGUMENT, &
         'handoff wrapper rejects a wrong trial extent')
    call check(short_trial(1) == 0.0_c_double .and. projected == 0_c_int &
         .and. handoff_is_zero(handoff_out), &
         'handoff extent failure clears trial, flag, and all ledger fields')
  end subroutine verify_bad_extents

end program test_fusion_handoff_fortran
