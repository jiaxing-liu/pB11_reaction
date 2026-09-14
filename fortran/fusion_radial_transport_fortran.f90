module fusion_radial_transport_fortran
  !! ISO_C_BINDING wrapper for the stateless radial transport trial.
  !!
  !! The first dimension of each rank-two array is contiguous in Fortran, so
  !! (zones, components), (zones+1, components), and (2, components) map
  !! directly to the component-major C buffers described by the API header.
  !! Every extent and the C-side entry bound is checked before C_LOC is used.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr, &
       c_size_t, c_loc
  implicit none
  private

  integer(c_int), parameter, public :: PB11_STATUS_OK = 0_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NULL_OUTPUT = 1_c_int
  integer(c_int), parameter, public :: PB11_STATUS_INVALID_ARGUMENT = 2_c_int
  integer(c_int), parameter, public :: PB11_STATUS_OUT_OF_RANGE = 3_c_int
  integer(c_int), parameter, public :: PB11_STATUS_NUMERICAL_FAILURE = 4_c_int
  integer(c_int), parameter, public :: PB11_STATUS_EXCEPTION = 5_c_int
  integer(c_int), parameter, public :: PB11_STATUS_UNKNOWN_METHOD = 6_c_int

  integer(c_size_t), parameter, public :: &
       FUSION_RADIAL_TRANSPORT_MAX_ENTRIES = 50000000_c_size_t

  ! Exact C layout: five consecutive c_double fields (40 bytes).
  type, bind(C), public :: fusion_radial_ledger_v1
     real(c_double) :: initial_number
     real(c_double) :: final_number
     real(c_double) :: inner_inward_number
     real(c_double) :: outer_outward_number
     real(c_double) :: balance_error
  end type fusion_radial_ledger_v1

  public :: fusion_radial_transport_trial

  interface
     function c_fusion_radial_transport_trial(zones, components, dt_s, &
          volume_old, volume_new, advection, conductance, old_density, &
          boundary_density, trial_density, face_amount, ledger) bind(C, &
          name="fusion_c_radial_transport_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: zones
       integer(c_int), value :: components
       real(c_double), value :: dt_s
       type(c_ptr), value :: volume_old
       type(c_ptr), value :: volume_new
       type(c_ptr), value :: advection
       type(c_ptr), value :: conductance
       type(c_ptr), value :: old_density
       type(c_ptr), value :: boundary_density
       type(c_ptr), value :: trial_density
       type(c_ptr), value :: face_amount
       type(c_ptr), value :: ledger
       integer(c_int) :: status
     end function c_fusion_radial_transport_trial
  end interface

contains

  subroutine clear_ledgers(ledger)
    type(fusion_radial_ledger_v1), intent(out) :: ledger(:)

    ledger%initial_number = 0.0_c_double
    ledger%final_number = 0.0_c_double
    ledger%inner_inward_number = 0.0_c_double
    ledger%outer_outward_number = 0.0_c_double
    ledger%balance_error = 0.0_c_double
  end subroutine clear_ledgers

  subroutine fusion_radial_transport_trial(zones, components, dt_s, &
       volume_old_m3, volume_new_m3, advection_volume_m3_s, &
       conductance_m3_s, old_density_m3, boundary_density_m3, &
       trial_density_m3, face_amount, ledger, status)
    integer(c_int), intent(in) :: zones, components
    real(c_double), intent(in) :: dt_s
    real(c_double), intent(in), target, contiguous :: volume_old_m3(:)
    real(c_double), intent(in), target, contiguous :: volume_new_m3(:)
    real(c_double), intent(in), target, contiguous :: advection_volume_m3_s(:)
    real(c_double), intent(in), target, contiguous :: conductance_m3_s(:,:)
    real(c_double), intent(in), target, contiguous :: old_density_m3(:,:)
    real(c_double), intent(in), target, contiguous :: boundary_density_m3(:,:)
    real(c_double), intent(out), target, contiguous :: trial_density_m3(:,:)
    real(c_double), intent(out), target, contiguous :: face_amount(:,:)
    type(fusion_radial_ledger_v1), intent(out), target, contiguous :: ledger(:)
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: z_count, component_count, face_count
    integer(c_size_t) :: cell_entries, face_entries, boundary_entries
    integer(c_size_t) :: total_factor, total_entries, max_size
    type(c_ptr) :: volume_old_ptr, volume_new_ptr, advection_ptr
    type(c_ptr) :: conductance_ptr, old_density_ptr, boundary_density_ptr
    type(c_ptr) :: trial_density_ptr, face_amount_ptr, ledger_ptr

    ! Make every output deterministic before validating dimensions.  A
    ! malformed extent returns before any C address is formed.
    trial_density_m3 = 0.0_c_double
    face_amount = 0.0_c_double
    call clear_ledgers(ledger)
    status = PB11_STATUS_INVALID_ARGUMENT

    volume_old_ptr = c_null_ptr
    volume_new_ptr = c_null_ptr
    advection_ptr = c_null_ptr
    conductance_ptr = c_null_ptr
    old_density_ptr = c_null_ptr
    boundary_density_ptr = c_null_ptr
    trial_density_ptr = c_null_ptr
    face_amount_ptr = c_null_ptr
    ledger_ptr = c_null_ptr

    if (zones < 1_c_int .or. components < 1_c_int) return

    z_count = int(zones, c_size_t)
    component_count = int(components, c_size_t)
    max_size = huge(0_c_size_t)

    ! Build all products in c_size_t and guard each addition/product before
    ! evaluating it.  The C implementation uses components*(2*zones+1).
    if (z_count > (max_size - 1_c_size_t) / 2_c_size_t) return
    total_factor = 2_c_size_t * z_count + 1_c_size_t
    if (component_count > max_size / total_factor) return
    total_entries = total_factor * component_count
    if (total_entries > FUSION_RADIAL_TRANSPORT_MAX_ENTRIES) return

    if (z_count >= max_size) return
    face_count = z_count + 1_c_size_t
    if (component_count > max_size / z_count) return
    cell_entries = z_count * component_count
    if (component_count > max_size / face_count) return
    face_entries = face_count * component_count
    if (component_count > max_size / 2_c_size_t) return
    boundary_entries = 2_c_size_t * component_count

    if (size(volume_old_m3, kind=c_size_t) /= z_count) return
    if (size(volume_new_m3, kind=c_size_t) /= z_count) return
    if (size(advection_volume_m3_s, kind=c_size_t) /= face_count) return

    if (size(conductance_m3_s, 1, kind=c_size_t) /= face_count) return
    if (size(conductance_m3_s, 2, kind=c_size_t) /= component_count) return
    if (size(conductance_m3_s, kind=c_size_t) /= face_entries) return

    if (size(old_density_m3, 1, kind=c_size_t) /= z_count) return
    if (size(old_density_m3, 2, kind=c_size_t) /= component_count) return
    if (size(old_density_m3, kind=c_size_t) /= cell_entries) return

    if (size(boundary_density_m3, 1, kind=c_size_t) /= 2_c_size_t) return
    if (size(boundary_density_m3, 2, kind=c_size_t) /= component_count) return
    if (size(boundary_density_m3, kind=c_size_t) /= boundary_entries) return

    if (size(trial_density_m3, 1, kind=c_size_t) /= z_count) return
    if (size(trial_density_m3, 2, kind=c_size_t) /= component_count) return
    if (size(trial_density_m3, kind=c_size_t) /= cell_entries) return

    if (size(face_amount, 1, kind=c_size_t) /= face_count) return
    if (size(face_amount, 2, kind=c_size_t) /= component_count) return
    if (size(face_amount, kind=c_size_t) /= face_entries) return
    if (size(ledger, kind=c_size_t) /= component_count) return

    ! All valid arrays have at least one element.  Their first dimensions are
    ! the contiguous dimensions, matching the C component-major layout.
    volume_old_ptr = c_loc(volume_old_m3(1))
    volume_new_ptr = c_loc(volume_new_m3(1))
    advection_ptr = c_loc(advection_volume_m3_s(1))
    conductance_ptr = c_loc(conductance_m3_s(1,1))
    old_density_ptr = c_loc(old_density_m3(1,1))
    boundary_density_ptr = c_loc(boundary_density_m3(1,1))
    trial_density_ptr = c_loc(trial_density_m3(1,1))
    face_amount_ptr = c_loc(face_amount(1,1))
    ledger_ptr = c_loc(ledger(1))

    status = c_fusion_radial_transport_trial(zones, components, dt_s, &
         volume_old_ptr, volume_new_ptr, advection_ptr, conductance_ptr, &
         old_density_ptr, boundary_density_ptr, trial_density_ptr, &
         face_amount_ptr, ledger_ptr)
    if (status /= PB11_STATUS_OK) then
       trial_density_m3 = 0.0_c_double
       face_amount = 0.0_c_double
       call clear_ledgers(ledger)
    end if
  end subroutine fusion_radial_transport_trial

end module fusion_radial_transport_fortran
