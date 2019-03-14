module time_mod
  use healpix_types
  use params
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values
  real(dp),    allocatable, dimension(:) :: H, z, omegaR, omegaLambda, omegaB, omegaM   ! Added grids

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point
  real(dp),                 dimension(1) :: y 
                                                                                                                                          

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2, n_total
    real(dp)     :: x_start, z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init, eps

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination
    n2          = 300                       ! Number of grid points after recombination
    n_t         = n1 + n2                   ! Number of grid points from beginning of recombination and onwards
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
   
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation
    eps         = 1.d-8                     ! error tolerance in ode_solver function


    
    ! Task: Fill in x and a grids
    allocate(x_t(n_t))
    allocate(a_t(n_t))
    allocate(H(n_t))
    allocate(z(n_t))
    allocate(omegaR(n_t))
    allocate(omegaLambda(n_t))
    allocate(omegaB(n_t))
    allocate(omegaM(n_t))

    dx = (x_end_rec - x_start_rec)/(n1 - 1)
    do i=1,n1
       x_t(i) = dx*(i - 1) + x_start_rec
    end do
    
    dx = (x_0 - x_end_rec)/n2 
    do i=n1,n_t
       x_t(i) = dx*(i - n1) + x_end_rec
    end do
    
    do i=1,n_t
       H(i) = get_H(x_t(i))*(3.09d19)
       z(i) = exp(-x_t(i)) - 1
       omegaR(i) = omega_r*H_zero**2*(exp(-4*x_t(i))/(get_H(x_t(i)))**2)
       omegaLambda(i) = omega_lambda*H_zero**2*(1/get_H(x_t(i))**2)
       omegaB(i) = omega_b*H_zero**2*(exp(-3*x_t(i))/(get_H(x_t(i)))**2)
       omegaM(i) =omega_m*H_zero**2*(exp(-3*x_t(i))/(get_H(x_t(i)))**2)
    end do

    do i=1,n_t
       a_t(i) = exp(x_t(i))
    end do

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
       
    dx = (x_eta2 - x_eta1)/(n_eta - 1)
    do i=1,n_eta
       x_eta(i) = x_eta1 + dx*(i-1)
    end do
 
    ! Solve ODE-system
    eta(1) = (c*a_init)/(H_zero*sqrt(omega_r)) ! Initialbetingelse
    y(1) = eta(1)
    do i=1,n_eta-1
     ! odeint takes arguments (y, x_start, x_end, error_tolerance, dx_start, dx_min, derivs_subroutine, method, output_subroutine)
     call odeint(y, x_eta(i), x_eta(i+1), eps, dx / 10.d0, dx / 10000.d0, derivs, bsstep, output)
     eta(i+1) = y(1)
    end do

    !do i=1,n_eta
    !  eta(i) = eta(i)/(3.09d22)
    !end do
    
    

  ! Output to file
  open(54,file='milestone1_part1.dat')
  do i=1,n_t
     write(54,'(7(E17.8))') x_t(i), z(i), H(i), omegaR(i), omegaLambda(i), omegaB(i), omegaM(i)
  end do
  close(54)

  open(54,file='milestone1_part2.dat')
  do i=1,n_eta
     write(54,'(2(E17.8))') x_eta(i), eta(i)/(3.09d22)
  end do
  close(54)

  write(*,*) "All data for Milestone 1 written to file"

  deallocate(z) 
  
  ! spline low resolution data
  allocate(eta2(n_eta))
  ! spline takes arguments (x, y, dy/dx|_1, dy/dx|_n, d^2y/dx^2),
  ! where the last argument is the output of the subroutine.
  ! Setting the two first derivatives to 1e30 corresponds to 
  ! choosing the "natural spline". 
  call spline(x_eta, eta, 1d30, 1d30, eta2)

  end subroutine initialize_time_mod

  subroutine derivs(x, y, dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    dydx(1) = c/get_H_p(x)
  end subroutine derivs

  subroutine output(x, y)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
  end subroutine output


  ! Task: Write a function that computes H at given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H
      
    get_H = H_zero*sqrt((omega_b + omega_m)*exp(-3*x) + omega_r*exp(-4*x) + omega_lambda)  

  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    
    get_H_p = H_zero*sqrt((omega_b + omega_m)*exp(-x) + omega_r*exp(-2*x) + omega_lambda*exp(2*x))

  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p

    get_dH_p = (H_zero/2)*(-(omega_b + omega_m)*exp(-x) - 2*omega_r*exp(-2*x) + omega_lambda*exp(x))/(sqrt((omega_b + omega_m)*exp(-x) + omega_r*exp(-2*x) + omega_lambda*exp(x)))

  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

    ! interpolate
    ! splint is the function doing the spline interpolation. 
    ! it takes as argument the three already computed arrays
    ! x, y and d^2y/dx^2 and then a single value of x that 
    ! you want to interpolate y to.  
    ! There exist a similar function splint_deriv, which 
    ! takes in the same arguments and returns the first 
    ! derivative of y at the interpolation point. 
    get_eta = splint(x_eta, eta, eta2, x_in)
    
  end function get_eta

end module time_mod
