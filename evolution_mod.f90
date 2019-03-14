module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
  implicit none

  ! Collecting terms
  real(dp) :: konstant2, konstant3, dx

  ! Accuracy parameters 
  real(dp),     parameter, private :: x_start_rec    = -log(1.d0 + 1630.4d0)   
  real(dp),     parameter, private :: a_init         = 1.d-8
  real(dp),     parameter, private :: x_init         = log(a_init)
  real(dp),     parameter, private :: k_min          = 0.1d0 * H_zero / c
  real(dp),     parameter, private :: k_max          = 1.d3  * H_zero / c
  integer(i4b), parameter          :: n_k            = 100
  integer(i4b), parameter, private :: lmax_int       = 6
    

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Storing TC-values
  real(dp), allocatable, dimension(:,:)   :: TC_verdier

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int
  
contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k, x
    real(dp), pointer, dimension(:,:), intent(out) :: S

    integer(i4b) :: i, j, n_hires
    real(dp)     :: g, dg, ddg, tau, dtau, ddtau, Hp, dHp, ddHHp, p, q, ddTheta2, x_tc, x_value, a, Theta0, Theta1, Theta2, Theta3, dTheta1, dTheta2, dTheta3, dx_hires, dk_hires, k_value
    real(dp), allocatable, dimension(:,:) :: S_lores

    real(dp), allocatable, dimension(:,:,:,:) :: coeff

    ! Task: Output a pre-computed 2D array (over k and x) for the                                                                                                                                                 
    !       source function, S(k,x). Remember to set up (and allocate) output                                                                                                                                     
    !       k and x arrays too.                                                                                                                                                                                   
    !                                                                                                                                                                                                             
    ! Substeps:                                                                                                                                                                                                   
    !   1) First compute the source function over the existing k and x                                                                                                                                                !      grids  

    allocate(S_lores(n_k,n_t))

    do j=1,n_k
       k_current = ks(j)
       do i=1,n_t
          x_value = x_t(i)
          a = exp(x_value)
          g = get_g(x_value)
          dg = get_dg(x_value)
          ddg = get_ddg(x_value)
          tau = get_tau(x_value)
          dtau = get_dtau(x_value)
          ddtau = get_ddtau(x_value)
          Hp = get_H_p(x_value)
          dHp = get_dH_p(x_value)
          ddHHp = (H_zero**2/2)*((omega_b + omega_m)*a**(-1) + 4*omega_r*a**(-2) + 4*omega_lambda*a**2)
          Theta0 = Theta(i,0,j)
          Theta1 = Theta(i,1,j)
          Theta2 = Theta(i,2,j)
          Theta3 = Theta(i,3,j)
          dTheta1 = dTheta(i,1,j)
          dTheta2 = dTheta(i,2,j)
          dTheta3 = dTheta(i,3,j)

          ddTheta2 = (2.d0*c*k_current/(5.d0*Hp))*(-dHp*Theta1/Hp + dTheta1) + (3.d0/10.d0)*(ddtau*Theta2 + dtau*dTheta2) - (3.d0*c*k_current/(5.d0*Hp))*(-dHp*Theta3/Hp + dTheta3)
          p = ddHHp*g*Theta2 + 3*Hp*dHp*(dg*Theta2 + g*dTheta2) + Hp**2*(ddg*Theta2 + 2*dg*dTheta2 + g*ddTheta2)       
          q = g*v_b(i,j)*dHp + Hp*v_b(i,j)*dg + Hp*g*dv_b(i,j)

          S_lores(j,i) = g*(Theta0 + Psi(i,j) + 0.25d0*Theta2) + exp(-tau)*(dPsi(i,j) - dPhi(i,j)) - (1.d0/(c*k_current))*q + (3.d0/(4.d0*c**2*k_current**2))*p
       end do
    end do

    !   2) Then spline this function with a 2D spline
    allocate(coeff(4,4,n_k,n_t))
    
    call splie2_full_precomp(ks,x_t,S_lores,coeff)

    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays
    
    n_hires = 5000
    allocate(x(n_hires))
    allocate(k(n_hires))
    allocate(S(n_hires,n_hires))
 
    dx_hires = (0 - x_start_rec)/(n_hires - 1)
    dk_hires = (k_max - k_min)/(n_hires - 1)

    do i=1,n_hires-1
       x(i) = x_start_rec + dx_hires*(i-1)
    end do
    x(n_hires) = 0.d0
    
    do j=1,n_hires
       k(j) = k_min + dk_hires*(j-1) 
    end do

    do j=1,n_hires
       k_value = k(j)
       do i=1,n_hires
          S(j,i) = splin2_full_precomp(ks,x_t,coeff,k_value,x(i))
       end do
    end do
    
  end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, n

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do i =1,n_k
       ks(i) = k_min + (k_max - k_min)*((i-1)/99.d0)**2
    end do
     
    n = 500
    ! Allocate array for storing TC-values
    allocate(TC_verdier(8,n))

    ! Allocate arrays for perturbation quantities
    allocate(Theta(0:n_t, 0:lmax_int, n_k))
    allocate(delta(0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(0:n_t, n_k))
    allocate(v_b(0:n_t, n_k))
    allocate(Phi(0:n_t, n_k))
    allocate(Psi(0:n_t, n_k))
    allocate(dPhi(0:n_t, n_k))
    allocate(dPsi(0:n_t, n_k))
    allocate(dv_b(0:n_t, n_k))
    allocate(dTheta(0:n_t, 0:lmax_int, n_k))

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)     = 1.d0
    delta(0,:)   = (3.d0/2.d0)
    delta_b(0,:) = (3.d0/2.d0)
    Theta(0,0,:) = 0.5d0
       
    do i = 1, n_k
       v(0,i)       = (c*ks(i))/(2*get_H_p(x_init))
       v_b(0,i)     = (c*ks(i))/(2*get_H_p(x_init))
       Theta(0,1,i) = -(c*ks(i))/(6*get_H_p(x_init))
       Theta(0,2,i) = -(20*c*ks(i)*Theta(0,1,i))/(45*get_H_p(x_init)*get_dtau(x_init))
       do l = 3, lmax_int
          Theta(0,l,i) = -(l/(2*l + 1.d0))*(c*ks(i)/(get_H_p(x_init)*get_dtau(x_init)))*Theta(0,l-1,i)
       end do
    end do

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l, n, n_x
    real(dp)     :: eps, x_tc, h, h1

    CHARACTER(LEN=128) :: number, filename

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx, x_grid, x_grid2
     
    h = 0.d0
    h1 = 1.d-5
    eps    = 1.d-8
    n = 500
   
    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))
    allocate(x_grid(0:n))

    ! Propagate each k-mode independently
    do k = 1, n_k

       k_current = ks(k)  ! Store k_current as a global module variable
       
       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k)
       y_tight_coupling(7) = Theta(0,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)
              
       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations

       ! New for each k
       konstant3 = (12*H_zero**2*omega_r)/(c**2*k_current**2)
       konstant2 = c*k_current

       ! Grid for integrating from x_init to x_tc
       dx = (x_tc - x_init)/n
       do i=0,n
          x_grid(i) = x_init + dx*i
       end do

       ! Solving set of differential equations during tight coupling. Final values will serve as the initial 
       ! values for solving differential equations from tight coupling up until the present epoch. 
       ! These initial conditions will be stored as the final values of y_tight_coupling
       do i=0,n-1
        
          call odeint(y_tight_coupling, x_grid(i), x_grid(i+1), eps, h1, h, derivs3, bsstep, output)
          do j=1,7
             TC_verdier(j,i+1) = y_tight_coupling(j)
          end do
          ! Psi
          TC_verdier(8,i+1) = -y_tight_coupling(5) + (konstant3/exp(2*x_grid(i+1)))*(20.d0/45.d0)*(konstant2/get_H_p(x_grid(i+1)))*(y_tight_coupling(7)/get_dtau(x_grid(i+1)))
       end do
        
       ! Writing TC-values to file, keeping file open
       if (k==1) then
          write(number,'(i1)') k
          filename = 'milestone3_k'//trim(number)//'.dat'
          open(54,file = filename)
          do i=1,n
             write(54,'(9(E19.8E3))') x_grid(i), TC_verdier(1,i), TC_verdier(2,i), TC_verdier(3,i), TC_verdier(4,i), TC_verdier(5,i),TC_verdier(6,i), TC_verdier(7,i), TC_verdier(8,i)
          end do  
       end if
 
       if (k==10 .or. k==30 .or. k==50 .or. k==80) then
          write(number,'(i2)') k
          filename = 'milestone3_k'//trim(number)//'.dat'
          open(54,file = filename)
          do i=1,n
             write(54,'(9(E19.8E3))') x_grid(i), TC_verdier(1,i), TC_verdier(2,i), TC_verdier(3,i), TC_verdier(4,i), TC_verdier(5,i),TC_verdier(6,i), TC_verdier(7,i), TC_verdier(8,i)
          end do
       end if
   
       if (k==100) then
          write(number,'(i3)') k
          filename = 'milestone3_k'//trim(number)//'.dat'
          open(54,file = filename)
          do i=1,n
             write(54,'(9(E19.8E3))') x_grid(i), TC_verdier(1,i), TC_verdier(2,i), TC_verdier(3,i), TC_verdier(4,i), TC_verdier(5,i),TC_verdier(6,i), TC_verdier(7,i), TC_verdier(8,i)
          end do
       end if 

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling(1:7) 
       y(8)   = -(20*c*k_current*y(7))/(45.d0*get_H_p(x_tc)*get_dtau(x_tc))
       do l = 3, lmax_int
          y(6+l) = -(l/(2*l + 1.d0))*(konstant2/(get_H_p(x_tc)*get_dtau(x_tc)))*y(6+l-1)
       end do
       
       n_x = ceiling(abs(x_start_rec - x_tc)/dx) 
       allocate(x_grid2(n_x))
       
       if (n_x > 1) then
          dx = (x_start_rec - x_tc)/(n_x - 1)
          do i=1,n_x
             x_grid2(i) = x_tc + dx*(i-1)
          end do
          do i=1,n_x-1
             call odeint(y, x_grid2(i), x_grid2(i+1), eps, h1, h, derivs4, bsstep, output)
          end do
       end if
       
       delta(1,k)   = y(1)
       delta_b(1,k) = y(2)
       v(1,k)       = y(3)
       v_b(1,k)     = y(4)
       Phi(1,k)     = y(5)
       do l = 0, lmax_int
          Theta(1,l,k) = y(6+l)
       end do
       Psi(1,k)     = - y(5) - konstant3*(y(8)/exp(2*x_t(1)))
                                                                                                                                           
       dPhi(1,k)     = dydx(5)
       dv_b(1,k)     = dydx(4)
       dTheta(1,:,k) = dydx(6:12)
       dPsi(1,k)     = - dydx(5) - (konstant3/exp(2*x_t(1)))*(dydx(8) - 2*y(8))

       deallocate(x_grid2)
                                                                                                                            
       do i = 1, n_t-1
          ! Task: Integrate equations from start of recombination to today
          call odeint(y, x_t(i), x_t(i+1), eps, h1, h, derivs4, bsstep, output)

          ! Task: Store variables at time step i in global variables
          delta(i+1,k)   = y(1)
          delta_b(i+1,k) = y(2)
          v(i+1,k)       = y(3)
          v_b(i+1,k)     = y(4)
          Phi(i+1,k)     = y(5)
          do l = 0, lmax_int
             Theta(i+1,l,k) = y(6+l)
          end do
          Psi(i+1,k)     = - y(5) - konstant3*(y(8)/exp(2*x_t(i+1)))

          ! Task: Store derivatives that are required for C_l estimation
          dPhi(i+1,k)     = dydx(5)
          dv_b(i+1,k)     = dydx(4)
          dTheta(i+1,:,k) = dydx(6:12)
          dPsi(i+1,k)     = - dydx(5) - (konstant3/exp(2*x_t(i+1)))*(dydx(8) - 2*y(8))
       end do
    
       if (k==1 .or. k==10 .or. k==30 .or. k==50 .or. k==80 .or. k==100) then
          do i=1,n_t
             write(54,'(9(E19.8E3))') x_t(i), delta(i,k), delta_b(i,k), v(i,k), v_b(i,k), Phi(i,k), Theta(i,0,k), Theta(i,1,k), Psi(i,k)
          end do
          close(54)
       end if
        
    end do
          
    write(*,*) "All data for Milestone 3 written to file"
    
    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)

  end subroutine integrate_perturbation_eqns

  subroutine derivs3(x,y,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp) :: a,R,andre_psi,q,theta2,Hp,dHp,dtau,ddtau,konstant1, phi_derivert, theta_null_derivert
    a = exp(x)
    Hp = get_H_p(x)
    dHp = get_dH_p(x)
    dtau = get_dtau(x)
    ddtau = get_ddtau(x)
    konstant1 = konstant2/Hp
    theta2 = -((20.d0/45.d0)*konstant1/dtau)*y(7)
    andre_psi = - y(5) - konstant3*(theta2/a**2)
    R = (4*omega_r)/(3*omega_b*a)
    
    phi_derivert = andre_psi - (konstant1**2/3)*y(5) + (H_zero**2/(2*Hp**2))*(omega_m*a**(-1.d0)*y(1) + omega_b*a**(-1.d0)*y(2) + 4*omega_r*a**(-2.d0)*y(6))
    theta_null_derivert = -konstant1*y(7) - phi_derivert
    q = (-((1-R)*dtau + (1+R)*ddtau)*(3*y(7) + y(4)) - konstant1*andre_psi + (1-dHp/Hp)*konstant1*(-y(6) + 2*theta2) - konstant1*theta_null_derivert) / ((1+R)*dtau + dHp/Hp - 1)

    dydx(1) = konstant1*y(3) - 3*phi_derivert
    dydx(2) = konstant1*y(4) - 3*phi_derivert
    dydx(3) = -y(3) - konstant1*andre_psi
    dydx(4) = (-y(4) - konstant1*andre_psi + R*(q + konstant1*(-y(6) + 2*theta2) - konstant1*andre_psi)) / (1+R)
    dydx(5) = phi_derivert
    dydx(6) = theta_null_derivert
    dydx(7) = (q - dydx(4))/3.d0
        
  end subroutine derivs3

  subroutine derivs4(x,y,dydx)
    use healpix_types
    implicit none
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx
    real(dp) :: a,R,andre_psi,Hp,dHp,dtau,ddtau,konstant1,phi_derivert
    a = exp(x)
    Hp = get_H_p(x)
    dHp= get_dH_p(x)
    dtau = get_dtau(x)
    ddtau = get_ddtau(x)
    konstant1 =konstant2/Hp
    R =(4*omega_r)/(3*omega_b*a)
    andre_psi = -y(5) - konstant3*(y(8)/a**2)

    phi_derivert = andre_psi - (konstant1**2/3)*y(5) + (H_zero**2/(2*Hp**2))*(omega_m*a**(-1.d0)*y(1) + omega_b*a**(-1.d0)*y(2) + 4*omega_r*a**(-2.d0)*y(6))    

    dydx(1) = konstant1*y(3) - 3*phi_derivert
    dydx(2) = konstant1*y(4) - 3*phi_derivert
    dydx(3) = -y(3) - konstant1*andre_psi
    dydx(4) = -y(4) - konstant1*andre_psi + dtau*R*(3*y(7) + y(4))
    dydx(5) = phi_derivert
    dydx(6) = -konstant1*y(7) - dydx(5) 
    dydx(7) = (konstant1/3.d0)*y(6) - (2.d0/3.d0)*konstant1*y(8) + (konstant1/3.d0)*andre_psi + dtau*(y(7) + y(4)/3.d0)
    dydx(8) = (2.d0/5.d0)*konstant1*y(7) - (3.d0/5.d0)*konstant1*y(9) + dtau*0.9d0*y(8)
    dydx(9) = (3.d0/7.d0)*konstant1*y(8) - (4.d0/7.d0)*konstant1*y(10) + dtau*y(9)
    dydx(10) = (4.d0/9.d0)*konstant1*y(9) - (5.d0/9.d0)*konstant1*y(11) + dtau*y(10)
    dydx(11) = (5.d0/11.d0)*konstant1*y(10) - (6.d0/11.d0)*konstant1*y(12) + dtau*y(11)
    dydx(12) = konstant1*y(11) - (7.d0*c*y(12))/(Hp*get_eta(x)) + dtau*y(12)
               
  end subroutine derivs4

  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none
    
    integer(i4b)          :: counter
    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time, factor1, factor2, x_value

    counter = 1
    factor1 = abs(get_dtau(x_init))
    factor2 = abs(c*k/(get_H_p(x_init)*factor1))
    x_value = x_init

    do while (factor1 > 10 .and. factor2 < 0.1d0 .and. x_value < x_start_rec)
       x_value = x_init - x_init*counter/10000.d0
       factor1 = abs(get_dtau(x_value))
       factor2 = abs(c*k/(get_H_p(x_value)*factor1))
       counter = counter + 1
    end do
    
    get_tight_coupling_time = x_value
       
  end function get_tight_coupling_time

end module evolution_mod
