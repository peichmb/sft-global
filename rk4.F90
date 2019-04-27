module rk4

    use params
    use grid
    use lsflows
    use vars
    use derivatives

    implicit none
    
    double precision, dimension(:,:), allocatable :: am1, am2, am3
    
contains

    subroutine init_rk4()
    
        integer :: i
    
        if (verbose) print *,'- Initializing rk4...'
    
        ! Allocate runge-kutta auxiliary matrices
        allocate(am1(nphi,nx)) ! x- advection  
        allocate(am2(nphi,nx)) ! phi- advection. ma
        
        ! Fill up auxiliary matrices
        do i = 1,nphi
            am1(i,:) = sqrt(1.d0-xx**2)/rsun ! SFT term 1
            am2(i,:) = -1.d0/rsun/sqrt(1.d0-xx**2) ! SFT term 2
        end do        
    
    end subroutine init_rk4

    subroutine pde_rhs(b, rhs, term)
    
    ! note the scope of b here!!
    ! this is not the same b as in the module vars; it is just used
    ! for the evaluation of the right hand side of the SFTE as needed
    ! by the Runge-Kutta method (b, b+k1/2 ...)
    
        double precision, dimension(nphi,nx) :: b, aux
        double precision, dimension(nphi,nx) :: bvtam1_x, bvp_phi, &
                                                & b_x
        double precision, dimension(nphi,nx) :: rhs

        integer :: i, term
        
        if (term == 1) then
            aux = b*vt*am1
            call deriv_x(aux,bvtam1_x)
            rhs = bvtam1_x
        else if (term == 2) then
            aux = b*vp
            call deriv_phi(aux,bvp_phi)
            rhs = am2*bvp_phi
        end if
            
    end subroutine pde_rhs
    
    subroutine rkstep(dtime,b,term)
    
        integer :: term
        double precision, intent(in) :: dtime
        double precision, dimension(nphi,nx) :: k1, k2, k3, k4, aux, b
               
        aux = b
        call pde_rhs(aux, k1, term)
        k1 = k1*dtime
        
        aux = b+k1/2.d0
        call pde_rhs(aux, k2, term)
        k2 = k2*dtime

        aux = b+k2/2.d0
        call pde_rhs(aux, k3, term)
        k3 = k3*dtime

        aux = b+k3
        call pde_rhs(aux, k4, term)
        k4 = k4*dtime
        
        b = b + (k1 + 2.d0*k2 + 2.d0*k3 + k4)/6.d0
        
    end subroutine rkstep
    
end module rk4
