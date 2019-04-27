module crank

    use params
    use grid
    use vars
    use timestep
    use derivatives
    use, intrinsic :: iso_c_binding
    
    implicit none
    
!    include 'fftw3.f'
    
    ! SFTE DIFFUSION TERM
    ! Variables for solution of 2nd order x-term
    double precision, dimension (:), allocatable :: alpha
    double complex, dimension(:), allocatable :: ldiag_x, udiag_x, diag0_x, diagm_x    
    ! Variables for solution of 2nd order phi-term
    double precision, dimension(:), allocatable :: beta, xi
    
    ! DIFFUSIVE SMOOTHING    
    ! Variables for solution of 2nd order x-term
    double precision, dimension (:), allocatable :: alpha_sm
    double complex, dimension(:), allocatable :: ldiag_x_sm, udiag_x_sm, diag0_x_sm, diagm_x_sm    
    ! Variables for solution of 2nd order phi-term
    double precision, dimension(:), allocatable :: beta_sm, xi_sm

contains
    
    subroutine init_crank()

        if (verbose) print *,'- Initializing crank...'
        
        ! DIFFUSION VARS
        ! --------------

        ! x- vars
        allocate(alpha(nx+1))
        allocate(ldiag_x(nx-1))
        allocate(diag0_x(nx))
        allocate(diagm_x(nx))
        allocate(udiag_x(nx-1))
        alpha = eta*dt0/rsun**2/dx**2*(1.d0-xx_12**2)/2.d0
        call diags_x(alpha, ldiag_x, diag0_x, diagm_x, udiag_x)
        
        ! phi- vars
        allocate(beta(nx))
        allocate(xi(nx))
        beta = eta*dt0/rsun**2/dphi**2/(1.d0-xx**2)
        xi = 2.d0*(1.d0+beta)
        
        
        ! DIFFUSIVE SMOOTHING
        ! -------------------
        
        ! x- vars
        allocate(alpha_sm(nx+1))
        allocate(ldiag_x_sm(nx-1))
        allocate(diag0_x_sm(nx))
        allocate(diagm_x_sm(nx))
        allocate(udiag_x_sm(nx-1))
        alpha_sm = eta*dt_sm/rsun**2/dx**2*(1.d0-xx_12**2)/2.d0
        call diags_x(alpha_sm, ldiag_x_sm, diag0_x_sm, diagm_x_sm, udiag_x_sm)

        ! phi- vars
        allocate(beta_sm(nx))
        allocate(xi_sm(nx))
        beta_sm = eta*dt_sm/rsun**2/dphi**2/(1.d0-xx**2)
        xi_sm = 2.d0*(1.d0+beta_sm)     

    end subroutine init_crank

!   ***********************************************************************    !
!   *****************                         *****************************    !
!   -----------------   Second order x-term   -----------------------------    !
!   *****************                         *****************************    !
!   ***********************************************************************    !

    subroutine diags_x(alpha, ldiag, diag0, diagm, udiag)

        double precision, dimension(nx+1) :: alpha
        double complex, dimension(nx-1) :: ldiag, udiag
        double complex, dimension (nx) :: diag0, diagm
        
        udiag = -alpha(2:nx)
        ldiag = -alpha(2:nx)
        diag0(2:nx-1) = 1.d0+alpha(2:nx-1)+alpha(3:nx)
        diagm(2:nx-1) = 1.d0+alpha(2:nx-1)+alpha(3:nx)

        ! Boundary conditions
        diag0(1) = 1.d0 + alpha(2)
        diag0(nx) = 1.d0 + alpha(nx)

        diagm(1) = 1.d0 + 2.d0*alpha(1) + alpha(2)
        diagm(nx) = 1.d0 + alpha(nx) + 2.d0*alpha(nx+1)
        
    end subroutine diags_x


    subroutine rhsvec_x(a, alpha, m, sol)
    
        integer :: m
        double precision, dimension(nx+1) :: alpha
        double complex, dimension(nx) :: a, sol

        sol(2:nx-1) = alpha(2:nx-1)*a(1:nx-2) &
                  & + (1.d0-alpha(2:nx-1)-alpha(3:nx))*a(2:nx-1) &
                  & + alpha(3:nx)*a(3:nx)
        
        if (m.eq.0) then
            sol(1) = (1.d0-alpha(2))*a(1) + alpha(2)*a(2)
            sol(nx) = alpha(nx)*a(nx-1) + (1.d0-alpha(nx))*a(nx)
        else
            sol(1) = (1.d0-2.d0*alpha(1)-alpha(2))*a(1) + alpha(2)*a(2)
            sol(nx) = alpha(nx)*a(nx-1) + (1.d0-alpha(nx)-2.d0*alpha(nx+1))*a(nx)
        end if
        
    end subroutine rhsvec_x

    subroutine idiffx(dtime, term)
    
        integer :: term
        double precision :: dtime
        
        if (term == 3) then
            alpha = alpha/dt_prev*dtime
            call diags_x(alpha, ldiag_x, diag0_x, diagm_x, udiag_x)
            call crankstepx(b, alpha, ldiag_x, diag0_x, diagm_x, udiag_x)
        else if (term == 101) then
            call crankstepx(b_sm, alpha_sm, ldiag_x_sm, diag0_x_sm, diagm_x_sm, udiag_x_sm)
        end if
        
    end subroutine idiffx

    
    subroutine crankstepx(b,alpha, ldiag, diag0, diagm, udiag)
        
        integer :: i, j, m, INFO, term
        integer*8 :: plan
        double complex, dimension(nx-1) :: udiag,ldiag, udiagc, ldiagc
        double complex, dimension(nx) :: diag0, diagm, diagc, sol
        double precision, dimension(nx+1) :: alpha
        double complex, dimension(nphi/2+1, nx) :: a
        double precision, dimension(nphi,nx) :: b
        
        ! Transform b to Fourier space
        do j = 1, nx
            call dfftw_plan_dft_r2c_1d(plan,nphi,b(:,j),a(:,j),FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c(plan,b(:,j),a(:,j))
            call dfftw_destroy_plan(plan)
        end do

        ! Solve tridiagonal systems
        do i = 1, nphi/2+1

            m = i-1
        
            ldiagc = ldiag
            udiagc = udiag
            if (m==0) then
                diagc = diag0
            else
                diagc = diagm
            end if
        
            call rhsvec_x(a(i,:), alpha, m, sol)
            
            call zgtsv(nx, 1, ldiagc, diagc, udiagc, sol, nx, INFO)
            
            a(i,:) = sol
            
        end do
                
        ! Transform back to real space
        do j = 1, nx
            call dfftw_plan_dft_c2r_1d(plan,nphi,a(:,j),b(:,j),FFTW_ESTIMATE)
            call dfftw_execute_dft_c2r(plan,a(:,j),b(:,j))
            call dfftw_destroy_plan(plan)
        end do
        
        b = b/float(nphi) ! FFTW scales the inverse transform
        
    end subroutine crankstepx



!   ***********************************************************************    !
!   *****************                         *****************************    !
!   -----------------  Second order phi-term  -----------------------------    !
!   *****************                         *****************************    !
!   ***********************************************************************    !
    
    
    subroutine diags_phi(beta, xi, ldiag, diag, udiag)

        double precision :: beta, xi
        double precision, dimension(nphi-1) :: ldiag, udiag
        double precision, dimension(nphi) :: diag

        udiag = -beta
        ldiag = -beta
        diag = xi
        diag(1) = 2.d0*xi
        diag(nphi) = xi+beta**2.d0/xi
                        
    end subroutine diags_phi


    subroutine rhsvec_phi(b,beta,sol)
        
        double precision :: beta
        double precision, dimension(nphi) :: b, sol
        
        sol(2:nphi-1) = beta*b(1:nphi-2) &
                    & + 2.d0*(1.d0-beta)*b(2:nphi-1) &
                    & + beta*b(3:nphi)
        
        ! Boundaries (periodic)
        sol(1) = beta*b(nphi) + 2.d0*(1.d0-beta)*b(1) + beta*b(2)
        sol(nphi) = beta*b(nphi-1) + 2.d0*(1.d0-beta)*b(nphi) + beta*b(1)
        
    end subroutine rhsvec_phi

    
    subroutine idiffp(dtime,term)
    
        double precision :: dtime
        integer :: term
        
        if(term == 4) then
            ! If dtime changes beta changes
            !if(inc_inflows) then
                beta = beta*dtime/dt_prev
                xi = 2.d0*(1.d0+beta)
                !print *, 1-xx(80)**2, dphi**2, rsun**2, eta, dtime, beta(80), xi(80)
            !end if
            call crankstepphi(b, beta, xi)
        else if(term == 102) then
            call crankstepphi(b_sm, beta_sm, xi_sm)            
        end if
    
    end subroutine idiffp
    
    
    subroutine crankstepphi(b, beta, xi)
        
        integer :: j, INFO
        double precision, dimension(nx) :: beta, xi
        double precision, dimension(nphi,nx) :: b
        double precision, dimension(nphi) :: sol, diag, y, z
        double precision, dimension(nphi-1) :: ldiag, udiag

        do j = 1, nx
            
            ! Solve auxiliary tridiagonal systems        
            call rhsvec_phi(b(:,j),beta(j),sol)           
            call diags_phi(beta(j), xi(j), ldiag, diag, udiag)
            y = sol
            call dgtsv(nphi, 1, ldiag, diag, udiag, y, nphi, INFO)

            call diags_phi(beta(j), xi(j), ldiag, diag, udiag)
            ! RHS vector is u = [-xi(j), 0, ..., 0, -beta(j)], see documentation.
            
            z = 0.d0
            z(1) = -xi(j)
            z(nphi) = -beta(j)
            call dgtsv(nphi, 1, ldiag, diag, udiag, z, nphi, INFO)
            
            ! Sherman-Morrison formula
            b(:,j) = y - ( y(1) + beta(j)/xi(j)*y(nphi) ) / &
                       & ( 1.d0 + z(1) + beta(j)/xi(j)*z(nphi) )*z
            
        end do

    end subroutine crankstepphi

    
    subroutine sediffp(b,eta,dtime)

        integer :: i, j
        integer*8 :: plan
        double precision, dimension(nphi,nx) :: b
        double complex, dimension(nphi/2+1) :: tf
        double precision :: dtime, eta


        do j = 1, nx
            
            ! Fourier transform
            call dfftw_plan_dft_r2c_1d(plan,nphi,b(:,j),tf,FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c(plan, b(:,j), tf)
            call dfftw_destroy_plan(plan)

            ! Diffuse in fourier space
            tf = tf*exp(-eta/sin2tht(j)/rsun**2*k_real**2*dtime)
            
            ! Transform back to real space
            call dfftw_plan_dft_c2r_1d(plan,nphi,tf,b(:,j),FFTW_ESTIMATE)
            call dfftw_execute_dft_c2r(plan, tf, b(:,j))
            call dfftw_destroy_plan(plan)
        
            b(:,j) = b(:,j)/float(nphi)
            
        end do

    end subroutine sediffp


end module crank
