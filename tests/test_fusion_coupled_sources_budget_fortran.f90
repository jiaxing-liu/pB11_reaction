program test_fusion_coupled_sources_budget_fortran
  use, intrinsic :: iso_c_binding
  use fusion_coupled_sources_budget_fortran
  implicit none
  integer(c_int), parameter :: cells=48_c_int, species=6_c_int, baths=7_c_int
  real(c_double), parameter :: ev_j=1.602176634e-19_c_double, mev_j=1.602176634e-13_c_double
  real(c_double), parameter :: dt_s=1.e-8_c_double
  real(c_double) :: edges(cells+1), thermal_number(species), thermal_charge_squared(species)
  real(c_double) :: coulomb_logs(baths,species), old_s(cells,species), old_t(cells,species)
  real(c_double) :: external_birth(cells,species), escape(cells,species)
  real(c_double) :: new_thermal_number(species), new_s(cells,species), new_t(cells,species)
  real(c_double) :: mapped(cells,7,5,2), electron_energy, ion_energy, electron_density
  type(fusion_coupled_thermal_options_v1) :: options
  type(fusion_fast_target_options_v1) :: fast_options
  type(fusion_coupled_thermal_v1), target :: output
  type(fusion_handoff_diagnostics_v1), target :: diagnostics
  type(fusion_beam_table_usage_v1), target :: usage
  type(fusion_coupled_floor_limits_v1) :: limits
  type(fusion_birth_floor_ledger_v1), target :: ledger
  type(fusion_birth_packets_v1), target :: packets
  type(fusion_coupled_source_accuracy_v1), target :: accuracy
  type(fusion_thermal_birth_budget_policy_v1) :: policy
  type(fusion_inert_ion_v1) :: inert(0)
  integer(c_int64_t) :: outside_count
  integer(c_int) :: status
  integer :: failures
  failures=0
  call check(c_sizeof(accuracy)==480_c_size_t, 'accuracy C layout is 480 bytes')
  call make_inputs()
  limits%max_center_over_ion_kT=1._c_double
  limits%max_ion_energy_fraction=1._c_double
  policy=fusion_thermal_birth_budget_policy_v1(1_c_int,4_c_int,4_c_int,0_c_int,0._c_double,0._c_double)
  call invoke(.false.)
  call check(status==PB11_STATUS_OK, 'valid inactive packet trial succeeds')
  call check(accuracy%version==1_c_int.and.accuracy%mode==0_c_int, 'inactive accuracy mode zero')
  call check(accuracy%number_L1_estimate==0._c_double.and.accuracy%energy_L1_estimate==0._c_double, &
       'inactive estimates are zero')
  call check(all(new_thermal_number==thermal_number).and.all(new_s==0._c_double), 'inactive state preserved')
  call seed()
  call invoke(.true.)
  call check(status==PB11_STATUS_INVALID_ARGUMENT, 'malformed packet shape rejected')
  call check_cleared()
  call seed()
  call invoke(.false.,.true.)
  call check(status==PB11_STATUS_INVALID_ARGUMENT, 'six-slot packet buffer rejected before C call')
  call check(all(mapped(:,:6,:,:)==0._c_double), 'six-slot rejected packet section cleared')
  call check(zero_bytes(c_loc(accuracy),c_sizeof(accuracy)), 'six-slot error clears accuracy')
  call check(all(new_thermal_number==0._c_double).and.outside_count==0_c_int64_t, &
       'six-slot error clears state and counter')
  call seed()
  policy%version=0_c_int
  call invoke(.false.)
  call check(status/=PB11_STATUS_OK, 'invalid policy rejected even with inactive pB')
  call check_cleared()
  call make_inputs()
  call make_pb_inputs()
  policy=fusion_thermal_birth_budget_policy_v1(1_c_int,4_c_int,4_c_int,0_c_int,2._c_double,2._c_double)
  call invoke(.false.)
  call check(status==PB11_STATUS_OK, 'live pB budget direct trial succeeds')
  call check(accuracy%version==1_c_int.and.accuracy%mode==1_c_int, 'live direct accuracy metadata')
  call check(accuracy%direct%version==1_c_int.and.accuracy%direct%mixed==1_c_int.and. &
       accuracy%direct%low_pairs>0_c_int64_t.and.accuracy%direct%high_pairs==0_c_int64_t, &
       'live report contains actual coarse parent counts')
  call check(accuracy%direct%full_nq==8_c_int.and.accuracy%direct%coarse_nq==4_c_int.and. &
       accuracy%direct%full_ncos==8_c_int.and.accuracy%direct%coarse_ncos==4_c_int, &
       'live report quadrature fields preserve C layout')
  call check(accuracy%number_L1_estimate>=accuracy%direct%number_L1_bound+ &
       accuracy%direct%number_roundoff_allowance.and.accuracy%energy_L1_estimate>= &
       accuracy%direct%energy_L1_bound+accuracy%direct%energy_roundoff_allowance, &
       'direct estimates include actual bounds and roundoff allowance')
  call check(output%ledger%events_m3(1)>0._c_double.and.sum(mapped(:,5,1,1))>0._c_double, &
       'live pB produces events and mapped alpha particles')
  call check(zero_bytes(c_loc(accuracy%table),c_sizeof(accuracy%table)), 'direct table info remains clear')
  if(failures/=0) error stop 1
  write(*,'(A)') 'Fortran coupled source budget wrapper checks passed'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition)then
      failures=failures+1
      write(*,'(A)') 'FAIL: '//message
    endif
  end subroutine
  subroutine invoke(bad_shape,short_species)
    logical,intent(in)::bad_shape
    logical,intent(in),optional::short_species
    integer :: last, slots
    slots=7
    if(present(short_species))then
      if(short_species)slots=6
    endif
    last=cells
    if(bad_shape)last=cells-1
    call fusion_coupled_sources_packets_budget_trial(dt_s,options,fast_options,cells,edges, &
         thermal_number,electron_energy,ion_energy,electron_density,thermal_charge_squared, &
         0_c_int,inert,coulomb_logs,old_s,old_t,external_birth,escape,new_thermal_number,new_s,new_t, &
         output,diagnostics,usage,limits,ledger,FUSION_BEAM_TABLE_STRICT,outside_count, &
         mapped(:last,:slots,:,:),packets,policy,accuracy,status)
  end subroutine
  subroutine seed()
    new_thermal_number=-7._c_double;new_s=-7._c_double;new_t=-7._c_double;mapped=-7._c_double
    outside_count=7_c_int64_t
    call fill_bytes(c_loc(output),c_sizeof(output),7_c_int8_t)
    call fill_bytes(c_loc(diagnostics),c_sizeof(diagnostics),7_c_int8_t)
    call fill_bytes(c_loc(usage),c_sizeof(usage),7_c_int8_t)
    call fill_bytes(c_loc(ledger),c_sizeof(ledger),7_c_int8_t)
    call fill_bytes(c_loc(packets),c_sizeof(packets),7_c_int8_t)
    call fill_bytes(c_loc(accuracy),c_sizeof(accuracy),7_c_int8_t)
  end subroutine
  subroutine fill_bytes(address,count,value)
    type(c_ptr),intent(in)::address
    integer(c_size_t),intent(in)::count
    integer(c_int8_t),intent(in)::value
    integer(c_int8_t),pointer::bytes(:)
    call c_f_pointer(address,bytes,[count])
    bytes=value
  end subroutine
  logical function zero_bytes(address,count)
    type(c_ptr),intent(in)::address
    integer(c_size_t),intent(in)::count
    integer(c_int8_t),pointer::bytes(:)
    call c_f_pointer(address,bytes,[count])
    zero_bytes=all(bytes==0_c_int8_t)
  end function
  subroutine check_cleared()
    call check(all(new_thermal_number==0._c_double).and.all(new_s==0._c_double).and. &
         all(new_t==0._c_double).and.outside_count==0_c_int64_t, 'failure clears state and count')
    ! The malformed-shape call owns only the section passed to the wrapper.
    call check(all(mapped(:cells-1,:,:,:)==0._c_double), 'failure clears mapped packet output')
    call check(zero_bytes(c_loc(output),c_sizeof(output)), 'failure clears complete coupled output')
    call check(zero_bytes(c_loc(diagnostics),c_sizeof(diagnostics)), 'failure clears diagnostics')
    call check(zero_bytes(c_loc(usage),c_sizeof(usage)), 'failure clears usage')
    call check(zero_bytes(c_loc(ledger),c_sizeof(ledger)), 'failure clears floor ledger')
    call check(zero_bytes(c_loc(packets),c_sizeof(packets)), 'failure clears packets')
    call check(zero_bytes(c_loc(accuracy),c_sizeof(accuracy)), 'failure clears every accuracy byte')
  end subroutine
  subroutine make_pb_inputs()
    integer :: j
    edges(1)=0._c_double
    do j=2,cells+1
      edges(j)=2.e-3_c_double*1.e3_c_double*ev_j* &
           (25._c_double*mev_j/(2.e-3_c_double*1.e3_c_double*ev_j))** &
           (real(j-2,c_double)/real(cells-1,c_double))
    enddo
    edges(cells+1)=25._c_double*mev_j
    thermal_number=0._c_double
    thermal_number(1)=5.e19_c_double
    thermal_number(6)=5.e19_c_double
    electron_density=3.e20_c_double
    electron_energy=1.5_c_double*electron_density*.5e3_c_double*ev_j
    ion_energy=1.5_c_double*sum(thermal_number)*.5e3_c_double*ev_j
    options%birth%relative_max_J=2._c_double*mev_j
    options%birth%cutoff_J=.001_c_double*mev_j
    options%birth%relative_order=4_c_int
    options%birth%cm_order=4_c_int
    options%max_source_rate_error=.05_c_double
    options%max_source_debit_error=.05_c_double
    options%channels(1)=1_c_int
    fast_options%angular_order=8_c_int
    limits%max_center_over_ion_kT=.5_c_double
    limits%max_ion_energy_fraction=.5_c_double
  end subroutine make_pb_inputs

  subroutine make_inputs()
    integer :: j

    do j = 1, cells + 1
       edges(j) = 25.0_c_double * mev_j * real(j - 1, c_double) / &
            real(cells, c_double)
    end do
    thermal_number = 0.0_c_double
    thermal_number(3) = 1.0e19_c_double
    thermal_charge_squared = [1.0_c_double, 1.0_c_double, 1.0_c_double, &
         4.0_c_double, 4.0_c_double, 25.0_c_double]
    coulomb_logs = 15.0_c_double
    old_s = 0.0_c_double
    old_t = 0.0_c_double
    external_birth = 0.0_c_double
    escape = 0.0_c_double
    electron_density = 1.0e19_c_double
    electron_energy = 1.5_c_double * electron_density * 5.0e3_c_double * ev_j
    ion_energy = 1.5_c_double * thermal_number(3) * 10.0e3_c_double * ev_j

    options%birth%relative_max_J = 5.0e6_c_double * ev_j
    options%birth%cm_max_kT = 40.0_c_double
    options%birth%ground_state_q_J = 91.84e3_c_double * ev_j
    options%birth%cutoff_J = 2.0e3_c_double * ev_j
    options%birth%l1_fraction = 0.76_c_double
    options%birth%relative_phase = 0.0_c_double
    options%birth%narrow_peak_fraction = 0.051_c_double
    options%birth%continuum_peak_scale = 1.0_c_double
    options%birth%continuation = 1_c_int
    options%birth%pb_low = 0_c_int
    options%birth%remainder_policy = 0_c_int
    options%birth%broad_mode = 13_c_int
    options%birth%fsci_policy = 0_c_int
    options%birth%relative_order = 8_c_int
    options%birth%cm_order = 8_c_int
    options%birth%nq = 8_c_int
    options%birth%ncos = 8_c_int
    options%max_source_rate_error = 2.0e-3_c_double
    options%max_source_debit_error = 2.0e-3_c_double
    options%handoff_max_L1 = 1.0e-3_c_double
    options%handoff_max_mean_error = 1.0e-3_c_double
    options%channels = 0_c_int
    options%handoff_enabled = 0_c_int

    fast_options%channels = 0_c_int
    fast_options%angular_order = 16_c_int
    fast_options%angular_max_exponent = 40.0_c_double
  end subroutine make_inputs

end program test_fusion_coupled_sources_budget_fortran
