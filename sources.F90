module sources

    use params
    use random
    use grid
    use vars
    use output
    
    implicit none
    
    integer :: icurr
    double precision :: dt_newar, rlatp, rlonp, rlatm, rlonm, rdb
    double precision, dimension(66) :: dbeta, dbetaprob ! size distribution
        
contains

    subroutine init_sources()
    
        integer :: i

        if (verbose) print *,'- Initializing sources...'
        
        ! New emergence probability is calculated for 1 day
        ! (i.e. number of new ARs per day)
        dt_newar = 1.d0*days2secs
        
        ! DISTRIBUTION FUNCTION N, SEE K. HARVEY, CHAPTER 12, (5)
        ! Number of AR with area A per day, n(A) ~ A^(-2)
        dbeta(1) = 3.5d0
        dbetaprob = 0.d0

        do i = 2, 66
    
            ! A ~ dbeta^2
            ! Size of ARs ranges from dbeta = 3.5 to dbeta = 10 deg
            dbeta(i) = 3.5d0 + float(i-1)/10.d0
       
            ! n(A) ~ A^(-2) --> n(A) ~ dbeta^(-4):
            ! Probability distr. must be normalized to unity (they sum up to 1)
            ! 14.04... makes dbetaprob(66) = 1
            dbetaprob(i) = dbetaprob(i-1) + 14.04628296399705d0/(dbeta(i)**4) 
            !dbetaprob(i) = dbetaprob(i-1) + 12.62d0/(dbeta(i)**4) ! IS THIS WRONG IN BAUMANN'S CODE?
       
        end do
        
        ! Active region record header
            
        open(123, file = trim(path)//'/ARrecord.'//suffix)
        write(123,'(A96)') '        Time    Latitude   Longitude    T. Angle     D. Beta        Area        Flux     U. flux'
        write(123,'(A96)') '      [days]       [deg]       [deg]       [deg]       [deg]     [deg^2]  [10^21 Mx]  [10^21 Mx]'
        write(123,*)
        close(123)

        if (armode ==2) then
            open(321, file = trim(inputrecord))
            read(321,*)
            read(321,*) icurr, rlatp, rlonp, rlatm, rlonm, rdb
        end if
            
    end subroutine init_sources

    
    function ncycles(time)
        
        double precision :: time, tau
        integer :: ncycles
        
        tau = time + tau0 - int((time+tau0)/cclmo)*cclmo
        
        if (tau < ccoverlap) then
            ncycles = 2
        else
            ncycles = 1
        end if
        
    end function ncycles
    
    
    function nnewars(tau)
    
        double precision :: tau, mean, x, y, actmod
        integer :: nnewars, n
        
        ! Gaussian distribution around mean number of new ARs
        mean = exp( -(tau-(cclength/2.d0))**2 / ((cclength/4.d0)**2) )
        
        x = rnd_u01()
        y = 1.26d0/mean
        
        n = 0
        if (x < 0.80d0/y) n = 4
        if (x < 0.78d0/y) n = 3
        if (x < 0.68d0/y) n = 2
        if (x < 0.48d0/y) n = 1

        ! Activity modifier
        actmod = n*actfac-int(n*actfac)
        n = int(n*actfac)
        x = rnd_u01()
        if (x < actmod) n = n+1    
        
        nnewars = n
    
    end function nnewars
    
    
    function tcycle(time, icyc)
    
        double precision :: time, tau, tcycle
        integer :: icyc
        
        tau = time + tau0 - int((time+tau0)/cclmo)*cclmo
        
        if (icyc == 2) tcycle = cclmo + tau
        if (icyc == 1) tcycle = tau
    
    end function tcycle
    
    
    function getarsize()
        
        double precision :: x, getarsize
        integer :: i
        
        x = rnd_u01()
        
        do i = 1, 66
      
           if (x <+ dbetaprob(i)) then
              getarsize = dbeta(i)
              exit
           end if
        end do

    end function getarsize
    
    
    function getarlon()
    
        double precision :: x, getarlon
        
        x = rnd_u01()
        
        getarlon = x*360.d0
        
    end function getarlon
    
    
    function getarlat(tau)
    
        double precision :: tau, getarlat, mu, sigma
    
        mu = alatu + (alatl-alatu)/cclength*tau
        sigma = alatscatu + (alatscatl-alatscatu)/cclength*tau
    
        getarlat = rnd_nor(mu,sigma)
        
        ! Half of the ARs appear in the southern hemisphere
        if (rnd_u01() < 0.5d0) getarlat = -getarlat
        
    end function getarlat
        
    
    function getartilt(lat)
    
        double precision :: getartilt, lat
        
        ! Gaussian distribution about the angle given by joy's law
        getartilt = rnd_nor(joysfac*lat,angscat)

    end function getartilt
    

    subroutine getarcoords(arlat, arlon, db, tilt, norm, thtp, phip, thtm, phim)
    
        logical :: norm
        double precision :: arlon, arlat, db, tilt, thtp, phip, thtm, phim
        double precision :: lonp, latp, lonm, latm

        arlat = arlat*deg2rad
        arlon = arlon*deg2rad
        tilt = tilt*deg2rad
        db = db*deg2rad
    
        if (arlat >= 0) then
            if (norm) then
                latp = arlat - db/2.d0*sin(tilt)
                lonp = arlon + db/2.d0*cos(tilt)/cos(arlat)
                latm = arlat + db/2.d0*sin(tilt)
                lonm = arlon - db/2.d0*cos(tilt)/cos(arlat)
            else
                latp = arlat + db/2.d0*sin(tilt)
                lonp = arlon - db/2.d0*cos(tilt)/cos(arlat)
                latm = arlat - db/2.d0*sin(tilt)
                lonm = arlon + db/2.d0*cos(tilt)/cos(arlat)
            end if
        else
            if (norm) then
                latp = arlat + db/2.*sin(tilt)
                lonp = arlon - db/2.*cos(tilt)/cos(arlat)
                latm = arlat - db/2.*sin(tilt)
                lonm = arlon + db/2.*cos(tilt)/cos(arlat)
            else
                latp = arlat - db/2.d0*sin(tilt)
                lonp = arlon + db/2.d0*cos(tilt)/cos(arlat)
                latm = arlat + db/2.d0*sin(tilt)
                lonm = arlon - db/2.d0*cos(tilt)/cos(arlat)        
            end if
        end if
    
        thtp = pi/2.d0-latp
        thtm = pi/2.d0-latm
 
        phip = lonp-2.d0*pi*int(lonp/2.d0/pi)
        phim = lonm-2.d0*pi*int(lonm/2.d0/pi)    
        
        arlat = arlat*rad2deg
        arlon = arlon*rad2deg
        tilt = tilt*rad2deg
        db = db*rad2deg
        
    end subroutine getarcoords
    

    subroutine getarvars(latp, lonp, latm, lonm, db, arlat, arlon, artilt)
    
        double precision :: latp, lonp, latm, lonm, db, arlat, arlon, artilt

        arlat = 0.5d0*(latp+latm)
        arlon = 0.5d0*(lonp+lonm)
        artilt = asin( (latp-latm)/db*(lonm-lonp)/abs(lonp-lonm) )*rad2deg

    end subroutine getarvars


    subroutine cyclepol(time, icycle, norm)

        logical :: norm
        double precision :: time
        integer :: icycle, cycleno
        
        norm = .true.
        cycleno = int((time+tau0)/cclmo)+1
        if (mod(cycleno,2) == 0) norm = .false.
        if (icycle ==2) norm = .not. norm
        if (.not. normstpol) norm = .not. norm
        
    end subroutine cyclepol
    
    subroutine mapar(thtp, phip, thtm, phim, db, bar)
    
        integer :: i, j
        double precision :: cosbetap, cosbetam, deltai2, delta02, edf
        double precision :: thtp, phip, thtm, phim, db, fluxp, fluxm
        double precision, dimension(nphi,ntht) :: barp, barm, bar

        db = db*deg2rad
        deltai2 = (0.4d0*db)**2
        delta02 = (4.d0*deg2rad)**2
        
        edf = deltai2/delta02 ! early diffusion factor; see Baumann (2004) 
            
        do j = 1, ntht
        do i = 1, nphi
            
            cosbetap = sin(thtp)*sintht(j)*cos(phip)*cosphi(i) &
                   & + sin(thtp)*sintht(j)*sin(phip)*sinphi(i) &
                   & + cos(thtp)*costht(j)
                
            barp(i,j) = bmax*edf*exp(-2.d0*(1.d0-cosbetap)/delta02)
                            
        end do
        end do
        
        do j = 1, ntht
        do i = 1, nphi
                   
            cosbetam = sin(thtm)*sintht(j)*cos(phim)*cosphi(i) &
                   & + sin(thtm)*sintht(j)*sin(phim)*sinphi(i) &
                   & + cos(thtm)*costht(j)
                
            barm(i,j) = -bmax*edf*exp(-2.d0*(1.d0-cosbetam)/delta02)
                            
        end do
        end do
        
        ! ARs appearing at high latitudes break flux conservation due to the poor
        !     resolution in theta near the poles. I need to slightly adjust the fluxes.
        fluxp = sum(barp)
        fluxm = sum(barm)        
        bar = barp-fluxp/fluxm*barm

        db = db*rad2deg
        
    end subroutine mapar


    subroutine placenewars(time)
    
        double precision, intent(in) :: time
        
        if (armode == 1) then
            call randomars(time)
        else if (armode == 2) then
            call filears(time)
        end if
        
    end subroutine placenewars

            
    subroutine filears(time)
    
        integer :: itime
        double precision, intent(in) :: time
        double precision :: thtp, phip, thtm, phim
        double precision :: arlat, arlon, artilt
        double precision, dimension(nphi,ntht) :: bar
        
        itime = nint(time*secs2days)        
        do while(icurr == itime)
            thtp = pi/2.d0 - rlatp*deg2rad
            phip = rlonp*deg2rad
            thtm = pi/2.d0 - rlatm*deg2rad
            phim = rlonm*deg2rad
            ! New AR map
            call mapar(thtp, phip, thtm, phim, rdb, bar)    
            ! Place new AR
            b = b + bar
            ! Add entry to record
            call getarvars(rlatp, rlonp, rlatm, rlonm, rdb, arlat, arlon, artilt)
            call addentry(time, arlat, arlon, artilt, rdb, bar)
            read(321,*) icurr, rlatp, rlonp, rlatm, rlonm, rdb
        end do

    end subroutine filears
    
    
    subroutine randomars(time)
    
        logical :: norm
        integer :: icyc, iar
        double precision, intent(in) :: time
        double precision :: tau, arlat, arlon, ardb, artilt
        double precision :: thtp, phip, thtm, phim, area, flux, uflux
        double precision, dimension(nphi,ntht) :: bar
        
        do icyc = 1, ncycles(time)
            
            tau = tcycle(time, icyc)

            do iar = 1, nnewars(tau)
                
                arlat = getarlat(tau)
                arlon = getarlon()
                ardb = getarsize()
                artilt = getartilt(arlat)
                
                call cyclepol(time, icyc, norm)
                call getarcoords(arlat, arlon, ardb, artilt, norm, thtp, phip, thtm, phim)
                
                ! New AR map
                call mapar(thtp, phip, thtm, phim, ardb, bar)
                
                ! Place new AR
                b = b + bar
                
                call addentry(time, arlat, arlon, artilt, ardb, bar)
                                                
            end do
            
        end do

    end subroutine randomars

    subroutine addentry(time, arlat, arlon, artilt, ardb, bar)

        double precision :: time, arlat, arlon, artilt, ardb, flux, uflux, area
        double precision, dimension(nphi,ntht) :: bar

        ! Calculate flux [10^21 Mx]
        flux = sum(bar)*(-dx)*dphi*rsuncgs**2/1.e21
        uflux = sum(abs(bar))*(-dx)*dphi*rsuncgs**2/1.e21
        
        ! Calculate area [deg^2]
        area = 2.d0*pi*(0.6d0*ardb)**2
                    
        ! Add entry to record
        open(123, position='append', file=trim(path)//'/ARrecord.'//suffix)
        write(123,'(8f12.3)') time*secs2days, arlat, arlon, artilt, ardb, area, flux, uflux
        close(123)

    end subroutine addentry

end module sources
