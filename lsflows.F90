module lsflows

    use params
    use grid
    
    implicit none
    
    double precision :: vm0
    double precision, dimension(:,:), allocatable :: mflow, drot
    
contains

    subroutine init_lsflows()
    
        if (verbose) print *,'- Initializing lsflows...'
        
        vm0 = 11.d0
    
        call fill_mflow()
        call fill_drot()
               
    end subroutine init_lsflows
    
    subroutine fill_mflow()
    
        integer :: j
        double precision :: lat
        
        allocate(mflow(nphi,nx))

        do j = 1, nx
            lat = pi/2.d0 - tt(j)
            if(abs(lat)*180.d0/pi > pcbnd) then
                mflow(:,j) = 0.d0
            else
                ! Note that lat is in rad but pcbnd is in deg.
                mflow(:,j) = - vm0*sin(180.*lat/pcbnd)
            end if
        end do
        
    end subroutine fill_mflow
    
    subroutine fill_drot()
    
        integer :: j
        double precision :: ddtrs, lat
        
        ddtrs = pi/180.d0/3600.d0/24.d0 ! deg/day to rad/s
        
        allocate(drot(nphi,nx))
        
        do j = 1, nx
            lat = pi/2.d0 - tt(j)
            drot(:,j) = (-2.30d0*sin(lat)**2 & 
                      & -1.62d0*sin(lat)**4 ) &
                      & *ddtrs & ! to rad/s
                      & *rsun*cos(lat) ! to m/s
        end do
        
    end subroutine fill_drot

end module lsflows
