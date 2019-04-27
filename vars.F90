module vars

    use params
    use grid
    
    implicit none
    
    double precision, dimension(:,:), allocatable :: b, vp, vt
    double precision, dimension(:,:), allocatable :: b_sm, vp_inf, vt_inf
    
contains

    subroutine init_vars()

        if (verbose) print *,'- Initializing vars...'
    
        ! SFT Equation
        ! ------------
        allocate(b(nphi,nx))
        allocate(vp(nphi,nx))
        allocate(vt(nphi,nx))
        
        b = 0.d0
        vp = 0.d0
        vt = 0.d0
        
        ! Inflows
        ! -------
        
        allocate(b_sm(nphi,nx))
        allocate(vt_inf(nphi,nx))
        allocate(vp_inf(nphi,nx))
        
        vp_inf = 0.d0
        vt_inf = 0.d0        
        
    end subroutine init_vars
    
end module vars
