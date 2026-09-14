program installed_table_coupling
 use, intrinsic :: iso_c_binding, only: c_int,c_double,c_ptr,c_null_ptr
 use fusion_birth_table_fortran, only: fusion_birth_table_control_v1, &
  fusion_birth_table_create,fusion_birth_table_destroy
 use fusion_coupled_thermal_fortran, only: fusion_coupled_thermal_options_v1, &
  fusion_inert_ion_v1,fusion_coupled_thermal_v1,fusion_coupled_thermal_table_trial
 implicit none
 integer(c_int),parameter::n=200_c_int
 real(c_double),parameter::keV=1.602176634e-16_c_double,MeV=1.602176634e-13_c_double,ne=1e22_c_double
 type(fusion_coupled_thermal_options_v1)::o
 type(fusion_birth_table_control_v1)::ctl
 type(fusion_coupled_thermal_v1)::out
 type(fusion_inert_ion_v1)::carbon(1)
 type(c_ptr)::handle,tables(5)
 real(c_double)::edges(n+1),zero(n,6),s(n,6),t(n,6),logs(8,6),Nion(6),nextN(6),z2(6),Ue,Ui
 integer(c_int)::status
 integer::j,low
 edges=0.0_c_double;zero=0.0_c_double;logs=15.0_c_double;low=n/3
 do j=1,low
  edges(j+1)=1e-10_c_double*keV*1e12_c_double**(real(j-1,c_double)/real(low-1,c_double))
 end do
 do j=low+1,n
  edges(j+1)=(100.0_c_double+24900.0_c_double*real(j-low,c_double)/real(n-low,c_double))*keV
 end do
 o%birth%relative_max_J=2.5_c_double*MeV;o%birth%cm_max_kT=40.0_c_double
 o%birth%ground_state_q_J=.09184_c_double*MeV;o%birth%cutoff_J=.001_c_double*MeV
 o%birth%l1_fraction=.76_c_double;o%birth%relative_phase=0.0_c_double
 o%birth%narrow_peak_fraction=.051_c_double;o%birth%continuum_peak_scale=1.0_c_double
 o%birth%continuation=1_c_int;o%birth%pb_low=0_c_int;o%birth%remainder_policy=0_c_int
 o%birth%broad_mode=13_c_int;o%birth%fsci_policy=0_c_int
 o%birth%relative_order=16_c_int;o%birth%cm_order=12_c_int;o%birth%nq=4_c_int;o%birth%ncos=4_c_int
 o%max_source_rate_error=1e-5_c_double;o%max_source_debit_error=1e-5_c_double
 o%handoff_max_L1=.001_c_double;o%handoff_max_mean_error=.001_c_double
 o%handoff_enabled=0_c_int;o%channels=0_c_int;o%channels(4)=1_c_int
 ctl%max_rate_error=.001_c_double;ctl%max_debit_error=.001_c_double
 ctl%max_number_L1=.001_c_double;ctl%max_energy_L1=.001_c_double
 ctl%max_direct_rate_discrepancy=1e-5_c_double;ctl%max_direct_debit_discrepancy=1e-5_c_double
 ctl%max_knots=64_c_int;ctl%max_evaluations=512_c_int;ctl%max_depth=12_c_int
 call fusion_birth_table_create(3_c_int,19.5_c_double*keV,23.0_c_double*keV,o%birth,ctl,edges,handle,status)
 if(status/=0)stop 1
 tables=c_null_ptr;tables(4)=handle
 carbon(1)%density_m3=.001_c_double*ne;carbon(1)%mass_kg=12.0_c_double*1.66053906892e-27_c_double
 carbon(1)%mean_charge_squared=36.0_c_double
 Nion=0.0_c_double;Nion(2)=(ne-6.0_c_double*carbon(1)%density_m3)/2.0_c_double;Nion(3)=Nion(2)
 z2=[1.0_c_double,1.0_c_double,1.0_c_double,4.0_c_double,4.0_c_double,25.0_c_double]
 Ue=1.5_c_double*ne*5.0_c_double*keV
 Ui=1.5_c_double*(sum(Nion)+carbon(1)%density_m3)*20.0_c_double*keV
 call fusion_coupled_thermal_table_trial(.000125_c_double,o,tables,n,edges,Nion,Ue,Ui,ne,z2, &
  1_c_int,carbon,logs,zero,zero,zero,zero,nextN,s,t,out,status)
 call fusion_birth_table_destroy(handle)
 if(status/=0)stop 2
 if(out%ledger%events_m3(4)<=0.0_c_double.or.out%inert_ion_heat_J_m3(5)<=0.0_c_double)stop 3
 if(nextN(2)>=Nion(2))stop 4
 print *, 'Installed Fortran table creation and actual DT coupled trial passed'
end program
