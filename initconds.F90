module initconds

    use params
    use grid
    use vars
    use lsflows
    use sources
    use inflows
    
    implicit none
    
contains

    subroutine set_initconds()

        if (verbose) print *,'- Setting initial conditions...'

        if (trim(scenario) == 'minnorm') then
        
            call minnorm()

        else if (trim(scenario) == 'zerofield') then
        
            call zerofield()

        else if (trim(scenario) == 'inflowtest') then
        
            call inflowtest()

        else if (trim(scenario) == 'inflowtest2') then
        
            call inflowtest2()

        else if (trim(scenario) == 'single') then
        
            call single()
            
        else if (trim(scenario) == 'dealias') then
        
            call dealias()
                                    
        else if (trim(scenario) == 'monopolar') then
        
            call monopolar()
                                    
        else
        
            print *, 'Initial condition not implemented'
            stop

        end if
    
    end subroutine set_initconds
    
    ! Minimum activity, normal polarity cycle
    subroutine minnorm()
        
        integer :: j
        double precision :: b0, a0, lat0
        
        !vp = drot
        !vt = mflow
        
        lat0 = pcbnd*deg2rad
        b0 = icondp1
        a0 = vm0*rsun*lat0/pi/eta

        do j = 1, ntht
            if(abs(lat(j)) >= lat0) then
                b(:,j) = lat(j)/abs(lat(j))*b0
            else
                b(:,j) = lat(j)/abs(lat(j))*b0*exp(-a0*(cos(pi*lat(j)/lat0)+1.d0))
            end if
        end do   
        
    end subroutine minnorm
    
    
    subroutine zerofield()
    
        b = 0.d0
        vp = drot
        vt = mflow
        
    end subroutine zerofield
    
    
    subroutine inflowtest()

        use crank
        use derivatives
        use output

        logical :: norm
        integer :: j,i
        double precision :: arlat, arlon, ardb, artilt
        double precision :: thtp, phip, thtm, phim
        double precision :: b0, a0, lat0
        double precision, dimension(nphi,ntht) :: bar
            
        vp = drot
        vt = mflow

        ! THIS IS NOT THE RIGHT PLACE TO INITIALIZE THESE MODULES
        ! I do it here just to be able to produce maps of the inflows
        ! but execution is intended to STOP at the end.
        ! If you remove the stop at the end of the routine you'll get a runtime
        ! severe error anyway.

        call init_crank()
        call init_derivatives()
        call init_output()

        arlat = icondp1
        arlon = 180.d0
        ardb = icondp2
        artilt = arlat/2.d0
        norm = .true.
        
        call getarcoords(arlat, arlon, ardb, artilt, norm, thtp, phip, thtm, phim)
        call mapar(thtp, phip, thtm, phim, ardb, bar)
        
        print *, thtp, thtm, phip, phim
        
        b = bar
        
        call add_inflows()
        print *, ' >>>> Max. inflow velocities (vt, vp, modulus): ', &
               & maxval(vt_inf), maxval(vp_inf), maxval(sqrt(vt_inf**2+vp_inf**2))
               
        call wrmaps(0.d0)
        
        stop
        
    end subroutine inflowtest
    
    subroutine inflowtest2()

        use crank
        use derivatives
        use output

        logical :: norm
        integer :: j,i
        double precision :: arlat, arlon, ardb, artilt, delta02
        double precision :: thtp, phip, thtm, phim
        double precision :: b0, a0, lat0, cosbetap
        double precision, dimension(nphi,ntht) :: bar
            
        vp = drot
        vt = mflow

        delta02 = (4.d0*deg2rad)**2
        ! THIS IS NOT THE RIGHT PLACE TO INITIALIZE THESE MODULES
        ! I do it here just to be able to produce maps of the inflows
        ! but execution is intended to STOP at the end.
        ! If you remove the stop at the end of the routine you'll get a runtime
        ! severe error anyway.

        call init_crank()
        call init_derivatives()
        call init_output()

        !thtp = pi/2.d0
        thtp = icondp2*deg2rad

        do j = 1, nx
            cosbetap = sin(thtp)*sintht(j)*cos(phip) &
                   & + sin(thtp)*sintht(j)*sin(phip) &
                   & + cos(thtp)*costht(j)
                   
            b(:,j) = bmax*exp(-2.d0*(1.d0-cosbetap)/delta02)
        end do

        call add_inflows()
        print *, ' >>>> Max. inflow velocities (vt, vp, modulus): ', &
               & maxval(vt_inf), maxval(vp_inf), maxval(sqrt(vt_inf**2+vp_inf**2))
               
        call wrmaps(0.d0)
        
        stop
        
    end subroutine inflowtest2


    subroutine single()
    
        logical :: norm
        integer :: j
        double precision :: arlat, arlon, ardb, artilt
        double precision :: thtp, phip, thtm, phim
        double precision, dimension(nphi,ntht) :: bar

        vp = drot
        vt = mflow
                
        arlat = icondp1
        arlon = 180.d0
        ardb = icondp2
        artilt = arlat/2.d0
        norm = .true.
        
        call getarcoords(arlat, arlon, ardb, artilt, norm, thtp, phip, thtm, phim)
        call mapar(thtp, phip, thtm, phim, ardb, bar)
        
        print *, sum(b)*(-dx)*dphi*rsuncgs**2/1.e21
        b = b + bar
        print *, sum(b)*(-dx)*dphi*rsuncgs**2/1.e21
        
        stop
        
    end subroutine single
    
    subroutine dealias()
    
        integer :: i, j
    
        b = 0.d0
        do j = 1, nx
        do i = 1, nphi
            if ( (nphi/2-i)**2+(nx/2-j)**2 <= 10) b(i,j) = 1000.d0
        end do
        end do
            
    end subroutine dealias
    
    subroutine monopolar()
    
        integer :: i, j
        double precision :: thtp, phip, cosbetap, delta02

        delta02 = (4.d0*deg2rad)**2
        phip = icondp1*deg2rad
        thtp = (90.d0-icondp2)*deg2rad
        
        b = 0.d0

        do j = 1, ntht
        do i = 1, nphi
            
            cosbetap = sin(thtp)*sintht(j)*cos(phip)*cosphi(i) &
                   & + sin(thtp)*sintht(j)*sin(phip)*sinphi(i) &
                   & + cos(thtp)*costht(j)
                   
            b(i,j) = bmax*exp(-2.d0*(1.d0-cosbetap)/delta02)
                            
        end do
        end do

            
    end subroutine monopolar    
    
end module initconds




