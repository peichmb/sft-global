module derivatives

    use params
    use grid
    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f'
    double precision, dimension(:), allocatable :: k, k_real, filter

contains


    subroutine init_derivatives()
    
        if (verbose) print *,'- Initializing derivatives...'
        
        allocate(k(nphi))
        allocate(k_real(nphi/2+1))
        allocate(filter(nphi))
        call wnums()
        
    
    end subroutine init_derivatives


    subroutine deriv_x(a, a_x)
    
        integer :: j
        integer*8 :: plan
        double precision, dimension(nphi,nx), intent(in) :: a
        double precision, dimension(nphi,nx), intent(out) :: a_x
        
        ! Differentiate inner cells
        do j = 3, nx-2
            a_x(:,j) = ( -a(:,j+2) + 8.d0*a(:,j+1)   &
                       & -8.d0*a(:,j-1) + a(:,j-2) ) &
                       & /12.d0/dx
        end do
        
        ! Differentiate cells near boundaries
        a_x(:,1) = ( -a(:,3) + 7.d0*a(:,2) + 8.d0*a(:,1) )/12.d0/dx
        a_x(:,2) = ( -a(:,4) + 8.d0*a(:,3) + -9.d0*a(:,1) )/12.d0/dx
        a_x(:,nx-1) = ( 9.d0*a(:,nx) - 8.d0*a(:,nx-2) + a(:,nx-3) )/12.d0/dx
        a_x(:,nx) = ( -8.d0*a(:,nx) - 7.d0*a(:,nx-1) + a(:,nx-2) )/12.d0/dx
        
    end subroutine deriv_x

subroutine deriv_x_fs(a, a_x)
    
        integer :: j
        integer*8 :: plan
        double precision, dimension(nphi,nx), intent(in) :: a
        double precision, dimension(nphi,nx), intent(out) :: a_x
        double precision, dimension(nphi) :: ift
        double complex, dimension(nphi/2+1) :: ft
        double complex, dimension(nphi/2+1, -1:nx+2) :: c
        double complex, dimension(nphi/2+1, nx) :: c_x
        
        ! Transform to Fourier space
        do j = 1, nx
            call dfftw_plan_dft_r2c_1d(plan,nphi,a(:,j),c(:,j),FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c(plan, a(:,j), c(:,j))
            call dfftw_destroy_plan(plan)
        end do
        
        ! Apply boundary conditions
        c(1,0) = c(1,1)
        c(1,-1) = c(1,2)
        c(2:,0) = -c(2:,1)
        c(2:,-1) = -c(2:,2)
        
        c(1,nx+1) = c(1,nx)
        c(1,nx+2) = c(1,nx-1)
        c(2:,nx+1) = -c(2:,nx)
        c(2:,nx+2) = -c(2:,nx-1)
        
        ! Differentiate
        do j = 1, nx
            c_x(:,j) = ( -c(:,j+2)+8.d0*c(:,j+1)- &
                       & 8.d0*c(:,j-1)+c(:,j-2) ) &
                       & /12.d0/dx
        end do
        
        ! Transform back to real space
        do j = 1, nx
            call dfftw_plan_dft_c2r_1d(plan,nphi,c_x(:,j),a_x(:,j),FFTW_ESTIMATE)
            call dfftw_execute_dft_c2r(plan,c_x(:,j),a_x(:,j))
            call dfftw_destroy_plan(plan)
        end do
        
        a_x = a_x/float(nphi)
        
    end subroutine deriv_x_fs

    
    subroutine deriv_phi(a, a_phi)
    
        double precision, dimension(nphi,nx), intent(in) :: a
        double precision, dimension(nphi,nx) :: a_phi
        double complex, dimension(nphi/2+1) :: tf
        integer :: i, j
        integer*8 :: plan

        do j = 1, nx
            
            ! Fourier transform
            call dfftw_plan_dft_r2c_1d(plan,nphi,a(:,j),tf,FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c(plan, a(:,j), tf)
            call dfftw_destroy_plan(plan)

            ! Derivative in fourier space
            tf = imunit*tf*k_real!*filter
            
            ! Transform back to real space
            call dfftw_plan_dft_c2r_1d(plan,nphi,tf,a_phi(:,j),FFTW_ESTIMATE)
            call dfftw_execute_dft_c2r(plan, tf, a_phi(:,j))
            call dfftw_destroy_plan(plan)
        
            a_phi(:,j) = a_phi(:,j)/float(nphi)
            
        end do
               
    end subroutine deriv_phi


    subroutine wnums()
    
        integer :: i
        
        k(1) = 0.d0
        k(nphi/2+1) = float(nphi)/2.d0
        do i = 2, nphi/2
            k(i) = float(i)-1
            k(nphi/2+i) = -float(nphi/2-i+1)
        end do
                
        k = 2.d0*pi*k/float(nphi)/dphi
        k_real = k(1:nphi/2+1)
        
        !filter = 1.d0
        !where(abs(k) >= 2.d0*maxval(k)/3.d0) filter = 0.d0
        !print *, k
        !print *, filter
        !stop
        
    end subroutine wnums
    

end module derivatives
