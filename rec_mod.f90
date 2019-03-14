module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n, n_hr                                      ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec, x_rec_hr                              ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22, tau_hr, dtau_hr, ddtau_hr  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2                                    ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22, g_hr, dg_hr, ddg_hr              ! Splined visibility function
  real(dp),                            private :: A, B, bigC, D, E, F                          ! Constants to make notation easier

contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i
    real(dp)     :: eps, saha_limit, dydx, xmin, xmax, dx, xstart, xstop, dx_hr 
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopp
    eps        = 1.d-8
    xmin       = 10.d0
    xmax       = 10000.d0
    n_hr       = 10000
    
    ! Collecting constants to ease notation
    A = (m_H/(omega_b*rho_c))*((m_e*k_b*T_0)/(2*pi*hbar*hbar))**(3.d0/2.d0)
    B = (0.448*64*pi*alpha**2*hbar*hbar*sqrt(epsilon_0)*(k_b*m_e*T_0)**(3.d0/2.d0))/(sqrt(27*pi*k_b*T_0)*m_e**2*c*(2*pi*hbar**2)**(3.d0/2.d0))
    bigC = epsilon_0/(k_b*T_0)
    D = -bigC/4.d0
    E = ((3*epsilon_0/(c*hbar))**3*m_H)/(64*pi**2*omega_b*rho_c)
    F = (0.448*64*pi*alpha**2*hbar**2*omega_b*rho_c*sqrt(epsilon_0))/(sqrt(27*pi*k_b*T_0)*m_e**2*c*m_H)

    allocate(x_rec(n))
    allocate(z(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    
    ! Task: Fill in x (rec) grid
    dx = (xstop - xstart)/(n-1)
    do i=1,n
       x_rec(i) = xstart + dx*(i-1)
       z(i) = exp(-x_rec(i)) - 1
    end do

    ! Task: Compute X_e and n_e at all grid times
    use_saha = .true.
    do i=1,n
       if (use_saha) then
          ! Use the Saha equation
          X_e(i) = 2.d0/(sqrt(1 + 4.d0/(A*get_f(x_rec(i)))) + 1)
          n_e(i) = log((omega_b*rho_c/m_H)*(X_e(i)/exp(3*x_rec(i))))
          if (X_e(i) < saha_limit) then
             use_saha = .false.
             y(1) = X_e(i)
          end if
       else
          ! Use the Peebles equation
          call odeint(y, x_rec(i-1), x_rec(i), eps, dx/xmin, dx/xmax, derivs1, bsstep, output)
          X_e(i) = y(1)
          n_e(i) = log((omega_b*rho_c/m_H)*(X_e(i)/exp(3*x_rec(i))))
       end if
    end do


    ! Task: Compute splined (log of) electron density function
    call spline(x_rec, n_e, 1d30, 1d30, n_e2)

    ! Task: Compute optical depth at all grid points
    tau(n) = exp(-20.d0)
    y(1) = tau(n)
    do i=1,n-1
       call odeint(y, x_rec(n+1-i), x_rec(n-i), eps, dx/xmin, dx/xmax, derivs2, bsstep, output)
       tau(n-i) = y(1) 
    end do

    do i=1,n
       tau(i) = log(tau(i)) ! Taking the log since it is easier to spline
    end do

    ! Task: Compute splined (log of) optical depth
    call spline(x_rec, tau, 1d30, 1d30, tau2)

    ! Task: Compute splined second derivative of (log of) optical depth
    call spline(x_rec, tau2, 1d30, 1d30, tau22)

    ! Task: Compute splined visibility function
    do i=1,n
       g(i) = -get_dtau(x_rec(i))*exp(-get_tau(x_rec(i)))
    end do

    call spline(x_rec, g, 1d30, 1d30, g2)

    ! Task: Compute splined second derivative of visibility function
    call spline(x_rec, g2, 1d30, 1d30, g22)

    ! Creating high-resolution arrays for plotting
    allocate(x_rec_hr(n_hr))
    allocate(tau_hr(n_hr))
    allocate(dtau_hr(n_hr))
    allocate(ddtau_hr(n_hr))
    allocate(g_hr(n_hr))
    allocate(dg_hr(n_hr))
    allocate(ddg_hr(n_hr))

    dx_hr = (xstop - xstart)/(n_hr - 1)

    do i=1,n_hr
       x_rec_hr(i) = xstart + dx_hr*(i-1)
       tau_hr(i) = get_tau(x_rec_hr(i))
       dtau_hr(i) = get_dtau(x_rec_hr(i))
       ddtau_hr(i) = get_ddtau(x_rec_hr(i))
       g_hr(i) = get_g(x_rec_hr(i))
       dg_hr(i) = get_dg(x_rec_hr(i))
       ddg_hr(i) = get_ddg(x_rec_hr(i))
    end do

    ! Writing to file
    open(54,file='milestone2_part1.dat')
    do i=1,n_hr
      write(54,'(7(E19.8E3))') x_rec_hr(i), tau_hr(i), dtau_hr(i), ddtau_hr(i), g_hr(i), dg_hr(i), ddg_hr(i)
    end do
    close(54)

    open(54,file='milestone2_part2.dat')
    do i=1,n
      write(54,'(2(E19.8E3))') z(i), X_e(i)
    end do
    close(54)


  write(*,*) "All data for Milestone 2 written to file"
 
  end subroutine initialize_rec_mod

  subroutine derivs1(x,y,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    dydx(1) = ((8.227d0 + get_H_p(x)*get_k(x)/(1-y(1)))/(8.227d0 + get_H_p(x)*get_k(x)/(1-y(1)) + get_small_h(x)))*((get_l(x)/get_H_p(x))*(1-y(1)) - (get_m(x)/get_H_p(x))*y(1)**2)
  end subroutine derivs1
 
  subroutine derivs2(x,y,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    dydx(1) = -(c*get_n_e(x)*sigma_T)/get_H(x)
  end subroutine derivs2

  

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    get_n_e = exp(splint(x_rec, n_e, n_e2, x))

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)    
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau
    get_tau = exp(splint(x_rec, tau, tau2, x))

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau
    get_dtau = exp(splint(x_rec, tau, tau2, x))*splint_deriv(x_rec, tau, tau2, x)

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    get_ddtau = exp(splint(x_rec, tau, tau2, x))*splint(x_rec, tau2, tau22, x) + (exp(splint(x_rec, tau, tau2, x))*splint_deriv(x_rec, tau, tau2, x))**2/exp(splint(x_rec, tau, tau2, x))

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g
    get_g = splint(x_rec, g, g2, x)

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg
    get_dg = splint_deriv(x_rec, g, g2, x)

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
    get_ddg = splint(x_rec, g2, g22, x)

  end function get_ddg

  function get_f(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_f
    get_f = exp(1.5d0*x)*exp(4*D*exp(x))

  end function get_f

  function get_k(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_k
    get_k = exp(2*x)*E

  end function get_k

  function get_small_h(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_small_h
    get_small_h = B*log(bigC*exp(x))*(exp(D*exp(x))/exp(x))

  end function get_small_h

  function get_l(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_l
    get_l = B*log(bigC*exp(x))*exp(4*D*exp(x))

  end function get_l

  function get_m(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_m
    get_m = (F*log(bigC*exp(x)))/exp(1.5d0*x)

  end function get_m

end module rec_mod
