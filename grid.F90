module grid

    use params
    
    implicit none
    
    ! Steps and coordinates
    double precision :: dphi, dx, dtht
    double precision, dimension(:), allocatable :: pp, tt, xx, lat, sin2tht
    double precision, dimension(:), allocatable :: xx_12
    double precision, dimension(:), allocatable :: sinphi, cosphi, sintht, costht
    double precision, dimension(:,:), allocatable :: st, ct, sl, cl, erfilter
    
contains

    subroutine init_grid()
    
        integer :: i, j

        if (verbose) print *,'- Initializing grid...'

        ! x and phi steps
        dphi = 2.d0*pi/float(nphi)
        dtht = pi/float(ntht)
        dx = -2.d0/float(nx)
        
        ! Allocate grid arrays
        ! REMINDER: nx = ntht 
        allocate(pp(nphi)) ! phi
        allocate(tt(ntht)) ! theta
        allocate(xx(nx)) ! x = cos(theta)
        allocate(xx_12(nx+1)) ! Points in between
        allocate(lat(nx)) ! latitude
        allocate(sinphi(nphi))
        allocate(cosphi(nphi))
        allocate(sintht(ntht))
        allocate(costht(ntht))
        allocate(sin2tht(ntht))
        allocate(st(nphi,nx))
        allocate(ct(nphi,nx))
        allocate(sl(nphi,nx))
        allocate(cl(nphi,nx))
        allocate(erfilter(nphi,nx))
        
        ! Fill up grid arrays
        do i = 1, nphi
            pp(i) = float(i-1)*dphi+dphi/2.d0
        end do

        do j = 1, nx
            xx(j) = 1.d0+float(j-1)*dx+dx/2.d0
            tt(j) = acos(xx(j))
            lat(j) = pi/2.d0-tt(j)
            st(:,j) = sin(tt(j))
            ct(:,j) = cos(tt(j))
        end do
        
        do i = 1, nphi
            erfilter(i,:) = 0.5d0*( erf((lat+latfilter*deg2rad)/wfilter/deg2rad) &
                                  - erf((lat-latfilter*deg2rad)/wfilter/deg2rad) )
        end do

        costht = xx
        sintht = sqrt(1.d0 - xx**2)        
        sin2tht = 1.d0 - xx**2
        cosphi = cos(pp)
        sinphi = sin(pp)
        xx_12(1:nx) = xx-dx/2.d0
        xx_12(nx+1) = xx_12(nx)+dx/2.d0
     
    end subroutine init_grid
     
end module grid

