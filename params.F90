module params

    implicit none
    
    ! Physical and mathematical constants
    double precision :: pi, rsun, rsuncgs
    double precision :: deg2rad, rad2deg, years2secs,days2secs, secs2years, secs2days
    double complex :: imunit
    
    ! Read-in parameters
    logical :: inc_inflows, wr_v, wr_b, verbose=.true., inc_newar, lfluxinf 
    logical :: wr_bfly, normstpol, latavinf, wr_adm, wr_cp
    integer :: nx, ntht, nphi, rseed, nsteps_sm, maxsteps, armode, lfiltermode
    character(100) :: path, scenario, inputrecord
    double precision :: icondp1, icondp2, icondp3
    double precision :: sigma_sm, dt_sm, fwhm_sm, latfilter, wfilter
    double precision :: pcbnd, eta, etakm2s, cfl, dt_write, t_final, dt_bfly, &
                     &  dt_avrg, dt_cpoint
    double precision :: cclength, ccoverlap, cclmo, joysfac, angscat, actfac
    double precision :: alatu, alatl, alatscatu, alatscatl, bmax
    double precision :: tau0, inf_a, inf_b
    
contains

    subroutine init_params(res)
        
        logical :: res
        
        if (verbose) print *,'- Initializing params...'
        
        ! Physical and mathematical constants      
        pi = 4.d0*atan(1.d0)
        deg2rad = pi/180.d0
        rad2deg = 180.d0/pi
        years2secs = 365.d0*24.d0*3600.d0
        secs2years = 1.d0/years2secs
        days2secs = 24.d0*3600.d0
        secs2days = 1.d0/days2secs
        rsun = 695.8d6 ! Meters
        rsuncgs = rsun*100.d0
        imunit = (0.d0,1.d0) ! Imaginary unit
        
        if (.not. res) call read_params()
        
        ! More parameters and conversions
        nx = ntht
        eta = etakm2s*1.e6
        t_final = t_final*years2secs
        dt_write = dt_write*days2secs
        dt_bfly = dt_bfly*days2secs
        cclength = cclength*years2secs
        ccoverlap = ccoverlap*years2secs
        cclmo = cclength - ccoverlap
        tau0 = tau0*years2secs
        dt_bfly = float(int(dt_bfly))
        dt_avrg = 1.d0*days2secs ! Save for average every 1 days
        dt_cpoint = dt_cpoint*years2secs   
        sigma_sm = rsun*2.d0*pi*fwhm_sm/360.d0/2.d0/sqrt(2.d0*log(2d0)) ! Meters
        dt_sm = sigma_sm**2/2.d0/eta/float(nsteps_sm)
        
    end subroutine init_params
    
    subroutine read_params()
    
        open(10,file='parameters')
        
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) ! --- Number of phi- and theta-points
        read(10,*) nphi
        read(10,*) ntht
        read(10,*) ! --- CFL parameter
        read(10,*) cfl
        read(10,*) ! --- Polar caps boundary (for meridional flow)
        read(10,*) pcbnd
        read(10,*) ! --- Horizontal diffusion coefficient [km^2/s]
        read(10,*) etakm2s
        read(10,*) ! --- Final time [years]
        read(10,*) t_final
        read(10,*) ! --- Maximum number of steps (0 = unlimited)
        read(10,*) maxsteps
        read(10,*) ! --- Random number generator seed
        read(10,*) rseed
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) ! --- Include inflows
        read(10,*) inc_inflows
        read(10,*) ! --- Azimuthally average inflows
        read(10,*) latavinf
        read(10,*) ! --- FWHM of smoothing gaussian [deg]
        read(10,*) fwhm_sm
        read(10,*) ! --- Number of substeps
        read(10,*) nsteps_sm
        read(10,*) ! --- a
        read(10,*) inf_a
        read(10,*) ! --- b
        read(10,*) inf_b
        read(10,*) ! --- Flux-based inflows
        read(10,*) lfluxinf
        read(10,*) ! --- Polar filter (0 = none, 1 = CJSS10, 2 = erf)
        read(10,*) lfiltermode
        read(10,*) ! --- Filter central latitude (+-)
        read(10,*) latfilter
        read(10,*) ! --- Filter width
        read(10,*) wfilter
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) ! --- Initial scenario
        read(10,*) scenario
        read(10,*) ! --- I. cond. parameter #1
        read(10,*) icondp1
        read(10,*) ! --- I. cond. parameter #2
        read(10,*) icondp2
        read(10,*) ! --- I. cond. parameter #3
        read(10,*) icondp3
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) ! --- Include new sources
        read(10,*) inc_newar
        read(10,*) ! --- AR input mode (1 = Random; 2 = Read from file)
        read(10,*) armode
        read(10,*) ! --- Input Active Region Record path
        read(10,*) inputrecord
        read(10,*) ! --- Cycle length [years]
        read(10,*) cclength
        read(10,*) ! --- Cycle overlap [years]
        read(10,*) ccoverlap
        read(10,*) ! --- Cycle time start [years]
        read(10,*) tau0
        read(10,*) ! --- Polarity of first cycle normal?
        read(10,*) normstpol
        read(10,*) ! --- Joy's law factor
        read(10,*) joysfac
        read(10,*) ! --- Tilt angle scattering
        read(10,*) angscat
        read(10,*) ! --- Activity factor
        read(10,*) actfac
        read(10,*) ! --- Mean active latitude at cycle start [degrees]
        read(10,*) alatu
        read(10,*) ! --- Latitude scattering at cycle start [degrees]
        read(10,*) alatscatu
        read(10,*) ! --- Mean active latitude at cycle end [degrees]
        read(10,*) alatl
        read(10,*) ! --- Latitude scattering at cycle end [degrees]
        read(10,*) alatscatl
        read(10,*) ! --- Maximum AR magnetic field density [gauss]
        read(10,*) bmax
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*) ! --- Verbose output
        read(10,*) verbose
        read(10,*) ! --- Output files path
        read(10,*) path
        read(10,*) ! --- Write velocities
        read(10,*) wr_v
        read(10,*) ! --- Write magnetic field
        read(10,*) wr_b
        read(10,*) ! --- Writing interval [days]
        read(10,*) dt_write
        read(10,*) ! --- Write azimuthally averaged b field
        read(10,*) wr_bfly
        read(10,*) ! --- Write axial dipole moment
        read(10,*) wr_adm
        read(10,*) ! -- Averaging / adm interval interval
        read(10,*) dt_bfly
        read(10,*) ! --- Write check points
        read(10,*) wr_cp
        read(10,*) ! -- Checkpointing interval
        read(10,*) dt_cpoint
        
        close(10)
    
    end subroutine read_params
    
    subroutine show_params()
    
        ! Show parameters
        if (verbose) then
            print *
            print *, '-- params.F90; show_params --'
            print *
            print *, '   Number of points (phi): ', nphi
            print *, '   Number of points (theta): ', ntht
            print *, '   CFL parameter: ', cfl
            print *, '   Polar caps boundary [deg]: ', pcbnd
            print *, '   Horizontal diffusion coefficient [km^2/s]: ', etakm2s
            print *, '   Final time [years]: ', t_final*secs2years
            print *, '   Random number generator seed: ',rseed
            print *
            print *, '   --- INFLOWS ---'
            print *, '   Include inflows?', inc_inflows
            if (inc_inflows) then
                print *, '   Azimuthally average inflows?', latavinf
                print *, '   FWHM of smoothing gaussian [deg]: ', fwhm_sm
                print *, '   Number of substeps: ', nsteps_sm
                print *, '   a: ', inf_a
                print *, '   b: ', inf_b
                print *, '   Light flux based inflows?: ', lfluxinf
                print *, '   Polar filter: ', lfiltermode
                print *, '   Filter central latitude [deg]: ', latfilter
                print *, '   Filter width [deg]: ', wfilter
            end if
            print *
            print *, '   --- INITIAL CONDITION ---'
            print *, '   Initial scenario: '//trim(scenario)
            print *, '   I. cond. parameter #1: ', icondp1
            print *, '   I. cond. parameter #2: ', icondp2
            print *, '   I. cond. parameter #3: ', icondp3
            print *
            print *, '   --- NEW SOURCES ---'
            print *, '   Include new sources? ', inc_newar
            if (inc_newar) then
                if (armode == 1) then
                    print *,'   Input mode: random.'
                    print *, '   Cycle length [years]: ', cclength*secs2years
                    print *, '   Cycle overlap [years]: ', ccoverlap*secs2years
                    print *, '   Cycle start time [years]: ', tau0*secs2years
                    print *, '   Is polarity of 1st cycle normal? ', normstpol
                    print *, "   Joy's law factor: ", joysfac
                    print *, '   Tilt angle scattering: ', angscat
                    print *, '   Activity factor: ', actfac
                    print *, '   Mean active latitude at cycle start [degrees]: ', alatu
                    print *, '   Latitude scattering at cycle start [degrees]: ', alatscatu
                    print *, '   Mean active latitude at cycle end [degrees]: ', alatl
                    print *, '   Latitude scattering at cycle end [degrees]: ', alatscatl
                else if (armode == 2) then
                    print *,'   Input mode: Read from '//trim(inputrecord)
                end if
                print *, '   Maximum AR magnetic field density [gauss]: ', bmax
            end if
            print *
            print *, '   --- OUTPUT ---'
            print *, '   Output interval [days]: ', dt_write*secs2days
            print *, '   Write inflow velocities? ', wr_v
            print *, '   Write magnetogram? ', wr_b
            print *, '   Output path: '//trim(path)
            print *, '   Verbose output? ', verbose
            print *, '   Write azimuthally averaged b field? ', wr_bfly
            print *, '   Write axial dipole moment?', wr_adm
            print *, '   Averaging interval ', dt_bfly*secs2days
            print *, '   Write check points? ', wr_cp
            print *, '   Checkpointing interval [years]', dt_cpoint*secs2years
            print *
        end if
    
    end subroutine show_params
        
end module params

