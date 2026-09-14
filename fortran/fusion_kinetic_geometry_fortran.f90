module fusion_kinetic_geometry_fortran
  !! ISO_C_BINDING wrapper for the stateless radial-plus-energy geometry trial.
  !!
  !! Densities are shaped (cells,6,zones), with the first dimension
  !! contiguous.  Conductances and boundaries use (face,cell,species) and
  !! (boundary,cell,species), respectively, so their flat Fortran layout is
  !! the component-major C layout used by fusion_c_kinetic_geometry_trial.
  !! Every extent and the C-side entry bound is checked before C_LOC is used.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr, &
       c_size_t, c_loc
  use fusion_source_volume_fortran, only : &
       fusion_transport_ledger_v1, PB11_STATUS_OK, PB11_STATUS_NULL_OUTPUT, &
       PB11_STATUS_INVALID_ARGUMENT, PB11_STATUS_OUT_OF_RANGE, &
       PB11_STATUS_NUMERICAL_FAILURE, PB11_STATUS_EXCEPTION, &
       PB11_STATUS_UNKNOWN_METHOD
  implicit none
  private

  integer(c_int), parameter, public :: &
       FUSION_KINETIC_GEOMETRY_MAX_ZONES = 300_c_int
  integer(c_int), parameter, public :: &
       FUSION_KINETIC_GEOMETRY_MAX_CELLS = 100000_c_int
  integer(c_size_t), parameter, public :: &
       FUSION_KINETIC_GEOMETRY_MAX_ENTRIES = 50000000_c_size_t

  public :: PB11_STATUS_OK
  public :: PB11_STATUS_NULL_OUTPUT
  public :: PB11_STATUS_INVALID_ARGUMENT
  public :: PB11_STATUS_OUT_OF_RANGE
  public :: PB11_STATUS_NUMERICAL_FAILURE
  public :: PB11_STATUS_EXCEPTION
  public :: PB11_STATUS_UNKNOWN_METHOD
  public :: fusion_transport_ledger_v1
  public :: fusion_kinetic_geometry_trial

  interface
     function c_fusion_kinetic_geometry_trial(zones, cells, dt_s, edges, &
          volume_old, volume_new, advection, compression, conductance, &
          boundary_s, boundary_t, old_s, old_t, trial_s, trial_t, &
          transport_step) bind(C, &
          name="fusion_c_kinetic_geometry_trial") result(status)
       import :: c_double, c_int, c_ptr
       integer(c_int), value :: zones
       integer(c_int), value :: cells
       real(c_double), value :: dt_s
       type(c_ptr), value :: edges
       type(c_ptr), value :: volume_old
       type(c_ptr), value :: volume_new
       type(c_ptr), value :: advection
       type(c_ptr), value :: compression
       type(c_ptr), value :: conductance
       type(c_ptr), value :: boundary_s
       type(c_ptr), value :: boundary_t
       type(c_ptr), value :: old_s
       type(c_ptr), value :: old_t
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       type(c_ptr), value :: transport_step
       integer(c_int) :: status
     end function c_fusion_kinetic_geometry_trial
  end interface

contains

  subroutine clear_transport_ledgers(ledger)
    type(fusion_transport_ledger_v1), intent(out) :: ledger(:)
    integer :: i

    do i = 1, size(ledger)
       ledger(i)%spatial_number = 0.0_c_double
       ledger(i)%spatial_energy_J = 0.0_c_double
       ledger(i)%work_J = 0.0_c_double
       ledger(i)%lower_number = 0.0_c_double
       ledger(i)%lower_energy_J = 0.0_c_double
       ledger(i)%upper_number = 0.0_c_double
       ledger(i)%upper_energy_J = 0.0_c_double
    end do
  end subroutine clear_transport_ledgers

  subroutine fusion_kinetic_geometry_trial(zones, cells, dt_s, edges_J, &
       volume_old_m3, volume_new_m3, advection_m3_s, compression_s_inv, &
       conductance_m3_s, boundary_s_m3, boundary_t_m3, old_s_m3, old_t_m3, &
       trial_s_m3, trial_t_m3, transport_step, status)
    integer(c_int), intent(in) :: zones, cells
    real(c_double), intent(in) :: dt_s
    real(c_double), intent(in), target, contiguous :: edges_J(:)
    real(c_double), intent(in), target, contiguous :: volume_old_m3(:)
    real(c_double), intent(in), target, contiguous :: volume_new_m3(:)
    real(c_double), intent(in), target, contiguous :: advection_m3_s(:)
    real(c_double), intent(in), target, contiguous :: compression_s_inv(:)
    real(c_double), intent(in), target, contiguous :: conductance_m3_s(:,:,:)
    real(c_double), intent(in), target, contiguous :: boundary_s_m3(:,:,:)
    real(c_double), intent(in), target, contiguous :: boundary_t_m3(:,:,:)
    real(c_double), intent(in), target, contiguous :: old_s_m3(:,:,:)
    real(c_double), intent(in), target, contiguous :: old_t_m3(:,:,:)
    real(c_double), intent(out), target, contiguous :: trial_s_m3(:,:,:)
    real(c_double), intent(out), target, contiguous :: trial_t_m3(:,:,:)
    type(fusion_transport_ledger_v1), intent(out), target, contiguous :: &
         transport_step(:)
    integer(c_int), intent(out) :: status

    integer(c_size_t) :: z_count, cell_count, species_count
    integer(c_size_t) :: face_count, component_count, density_entries
    integer(c_size_t) :: conductance_entries, boundary_entries
    integer(c_size_t) :: total_factor, total_entries, max_size
    type(c_ptr) :: edges_ptr, volume_old_ptr, volume_new_ptr
    type(c_ptr) :: advection_ptr, compression_ptr, conductance_ptr
    type(c_ptr) :: boundary_s_ptr, boundary_t_ptr, old_s_ptr, old_t_ptr
    type(c_ptr) :: trial_s_ptr, trial_t_ptr, transport_ptr

    ! Clear every output before checking dimensions.  Invalid extents return
    ! before any C address is formed, including zero-size actual arrays.
    trial_s_m3 = 0.0_c_double
    trial_t_m3 = 0.0_c_double
    call clear_transport_ledgers(transport_step)
    status = PB11_STATUS_INVALID_ARGUMENT

    edges_ptr = c_null_ptr
    volume_old_ptr = c_null_ptr
    volume_new_ptr = c_null_ptr
    advection_ptr = c_null_ptr
    compression_ptr = c_null_ptr
    conductance_ptr = c_null_ptr
    boundary_s_ptr = c_null_ptr
    boundary_t_ptr = c_null_ptr
    old_s_ptr = c_null_ptr
    old_t_ptr = c_null_ptr
    trial_s_ptr = c_null_ptr
    trial_t_ptr = c_null_ptr
    transport_ptr = c_null_ptr

    if (zones < 1_c_int .or. zones > FUSION_KINETIC_GEOMETRY_MAX_ZONES) return
    if (cells < 1_c_int .or. cells > FUSION_KINETIC_GEOMETRY_MAX_CELLS) return

    z_count = int(zones, c_size_t)
    cell_count = int(cells, c_size_t)
    species_count = 6_c_size_t
    max_size = huge(0_c_size_t)

    ! Guard each c_size_t product before evaluating it.  The C entry bound is
    ! 12*cells*(2*zones+1), i.e. two six-species components plus radial data.
    if (z_count > (max_size - 1_c_size_t) / 2_c_size_t) return
    total_factor = 2_c_size_t * z_count + 1_c_size_t
    if (cell_count > max_size / species_count) return
    component_count = species_count * cell_count
    if (component_count > max_size / 2_c_size_t) return
    if (2_c_size_t * component_count > max_size / total_factor) return
    total_entries = (2_c_size_t * component_count) * total_factor
    if (total_entries > FUSION_KINETIC_GEOMETRY_MAX_ENTRIES) return

    if (z_count >= max_size) return
    face_count = z_count + 1_c_size_t
    if (component_count > max_size / z_count) return
    density_entries = component_count * z_count
    if (component_count > max_size / face_count) return
    conductance_entries = component_count * face_count
    if (component_count > max_size / 2_c_size_t) return
    boundary_entries = 2_c_size_t * component_count

    if (size(edges_J, kind=c_size_t) /= cell_count + 1_c_size_t) return
    if (size(volume_old_m3, kind=c_size_t) /= z_count) return
    if (size(volume_new_m3, kind=c_size_t) /= z_count) return
    if (size(advection_m3_s, kind=c_size_t) /= face_count) return
    if (size(compression_s_inv, kind=c_size_t) /= z_count) return

    if (size(conductance_m3_s, 1, kind=c_size_t) /= face_count) return
    if (size(conductance_m3_s, 2, kind=c_size_t) /= cell_count) return
    if (size(conductance_m3_s, 3, kind=c_size_t) /= species_count) return
    if (size(conductance_m3_s, kind=c_size_t) /= conductance_entries) return

    if (size(boundary_s_m3, 1, kind=c_size_t) /= 2_c_size_t) return
    if (size(boundary_s_m3, 2, kind=c_size_t) /= cell_count) return
    if (size(boundary_s_m3, 3, kind=c_size_t) /= species_count) return
    if (size(boundary_s_m3, kind=c_size_t) /= boundary_entries) return
    if (size(boundary_t_m3, 1, kind=c_size_t) /= 2_c_size_t) return
    if (size(boundary_t_m3, 2, kind=c_size_t) /= cell_count) return
    if (size(boundary_t_m3, 3, kind=c_size_t) /= species_count) return
    if (size(boundary_t_m3, kind=c_size_t) /= boundary_entries) return

    if (size(old_s_m3, 1, kind=c_size_t) /= cell_count) return
    if (size(old_s_m3, 2, kind=c_size_t) /= species_count) return
    if (size(old_s_m3, 3, kind=c_size_t) /= z_count) return
    if (size(old_s_m3, kind=c_size_t) /= density_entries) return
    if (size(old_t_m3, 1, kind=c_size_t) /= cell_count) return
    if (size(old_t_m3, 2, kind=c_size_t) /= species_count) return
    if (size(old_t_m3, 3, kind=c_size_t) /= z_count) return
    if (size(old_t_m3, kind=c_size_t) /= density_entries) return
    if (size(trial_s_m3, 1, kind=c_size_t) /= cell_count) return
    if (size(trial_s_m3, 2, kind=c_size_t) /= species_count) return
    if (size(trial_s_m3, 3, kind=c_size_t) /= z_count) return
    if (size(trial_s_m3, kind=c_size_t) /= density_entries) return
    if (size(trial_t_m3, 1, kind=c_size_t) /= cell_count) return
    if (size(trial_t_m3, 2, kind=c_size_t) /= species_count) return
    if (size(trial_t_m3, 3, kind=c_size_t) /= z_count) return
    if (size(trial_t_m3, kind=c_size_t) /= density_entries) return
    if (size(transport_step, kind=c_size_t) /= z_count) return

    ! Every valid array has at least one element; only now form C addresses.
    edges_ptr = c_loc(edges_J(1))
    volume_old_ptr = c_loc(volume_old_m3(1))
    volume_new_ptr = c_loc(volume_new_m3(1))
    advection_ptr = c_loc(advection_m3_s(1))
    compression_ptr = c_loc(compression_s_inv(1))
    conductance_ptr = c_loc(conductance_m3_s(1,1,1))
    boundary_s_ptr = c_loc(boundary_s_m3(1,1,1))
    boundary_t_ptr = c_loc(boundary_t_m3(1,1,1))
    old_s_ptr = c_loc(old_s_m3(1,1,1))
    old_t_ptr = c_loc(old_t_m3(1,1,1))
    trial_s_ptr = c_loc(trial_s_m3(1,1,1))
    trial_t_ptr = c_loc(trial_t_m3(1,1,1))
    transport_ptr = c_loc(transport_step(1))

    status = c_fusion_kinetic_geometry_trial(zones, cells, dt_s, edges_ptr, &
         volume_old_ptr, volume_new_ptr, advection_ptr, compression_ptr, &
         conductance_ptr, boundary_s_ptr, boundary_t_ptr, old_s_ptr, &
         old_t_ptr, trial_s_ptr, trial_t_ptr, transport_ptr)
    if (status /= PB11_STATUS_OK) then
       trial_s_m3 = 0.0_c_double
       trial_t_m3 = 0.0_c_double
       call clear_transport_ledgers(transport_step)
    end if
  end subroutine fusion_kinetic_geometry_trial

end module fusion_kinetic_geometry_fortran
