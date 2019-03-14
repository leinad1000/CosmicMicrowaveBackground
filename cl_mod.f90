module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline, n_hires
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e, dz, x_start_rec, sum, k_value, dk, k_min, k_max, pi
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand, etas
    real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp
    real(dp),     pointer,     dimension(:)       :: k, x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     pointer,     dimension(:,:)     :: S, S2, S_hires
    real(dp),     allocatable, dimension(:,:)     :: Theta
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2, cls_hires
    real(dp),     pointer,     dimension(:)       :: x_hires, k_hires

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! Set up which l's to compute
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    ! Task: Get source function from evolution_mod
    n_hires = 5000
    !allocate(k_hires(n_hires))
    !allocate(x_hires(n_hires))
    !allocate(S_hires(n_hires,n_hires))
    call get_hires_source_function(k_hires, x_hires, S_hires)
    
    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))
    
    ! Array for dummy variable, z, of Bessel-function
    dz = 3500.d0/(n_spline - 1)
    do i=1,n_spline
       z_spline(i) = dz*(i-1)
    end do

 
    ! Bessel-function for all 44 l's over array of z for each l
    do l=1,l_num
       j_l(1,l) = 0
       do i=2,n_spline
          call sphbes(ls(l), z_spline(i), j_l(i,l))
       end do
    end do

    allocate(Theta(l_num,n_hires))
    allocate(cls(l_num))
    allocate(cls2(l_num))
    allocate(etas(n_hires))

    eta0 = 3.4d0*c/H_zero 
    x_start_rec = -log(1.d0 + 1630.4d0)
    dx = (0 - x_start_rec)/(n_hires - 1)
    k_min = 0.1d0*H_zero/c
    k_max = 1.d3*H_zero/c
    dk = (k_max - k_min)/(n_hires - 1)
    sum = 0.d0
    pi = 4*atan(1.d0)
    
    do i=1,n_hires
       etas(i) = eta0 - get_eta(x_hires(i))
    end do
    
    ! Overall task: Compute the C_l's for each given l
    do l = 1, l_num

       ! Task: Compute the transfer function
       call spline(z_spline,j_l(:,l),1d30,1d30,j_l2(:,l))
              
       do j=1,n_hires
          k_value = k_hires(j)
          do i=1,n_hires
             z = k_value*etas(i)
             S_func = S_hires(j,i)
             j_func = splint(z_spline,j_l(:,l), j_l2(:,l), z)
             sum = sum + S_func*j_func
          end do
          Theta(l,j) = sum*dx
          sum = 0.d0
       end do
       
       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
       do j=1,n_hires
          sum = sum + (c*k_hires(j)/H_zero)**(-0.04d0)*(Theta(l,j)**2/(k_hires(j)))
       end do
       
       ! Task: Store C_l in an array. Optionally output to file
       cls(l) = sum*dk*(ls(l)*(ls(l)+1.d0)/(2*pi))
       sum = 0.d0
    end do

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    call spline(real(ls,dp),cls,1d30,1d30,cls2)
    
    allocate(cls_hires(1200))

    do l=2,1200
       cls_hires(l) = splint(real(ls,dp),cls,cls2,real(l,dp))
    end do

    cls_hires(1) = cls_hires(2)
    
    ! Writing to file
    open(54,file='milestone4.dat')
    do l=1,1200
       write(54,'(2(E19.8E3))') real(l,dp), 146715.10593973883*cls_hires(l)
    end do
    close(54)
     
    write (*,*) "All data for Milestone 4 written to file"    
 
  end subroutine compute_cls
  
end module cl_mod
