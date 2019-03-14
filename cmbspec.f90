program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use cl_mod
  implicit none

  ! Initialize time grids
  call initialize_time_mod
  call initialize_rec_mod
  call initialize_perturbation_eqns
  call integrate_perturbation_eqns  
  call compute_cls

end program cmbspec
