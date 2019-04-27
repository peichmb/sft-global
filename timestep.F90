module timestep

    use params
    use vars
    use grid
    use sources
    
    implicit none
    
    double precision :: prdtime_maps, prdtime_bfly, prdtime_nar, prdtime_avrg, &
                        dt0, prdtime_cpoint, dt_prev
    double precision :: tacc_nar = 0.d0, tacc_bfly, tacc_maps = 0.d0, &
                      & tacc_avrg = 0.d0, tacc_cpoint = 0.d0
    
contains 

    subroutine init_timestep()
    
        call dtmax(dt0)
        dt_prev = dt0

        if (verbose) then
            print *,'- Initializing timestep...'
            print *
            print *, '-- timestep.F90; init_timestep() --'
            print *
            print *, '   dt0 = ',dt0,' s'
            print *
        end if
    
    end subroutine init_timestep

    subroutine calc_dt(dtime, twmaps, newars, twbfly, twavrg, twcpoint, time, t_final, endloop)
    
        logical :: twmaps, endloop, newars, twbfly, twavrg, twcpoint
        double precision :: dtime, time, t_final

        ! If inflows are included, the velocity is time-dependent,
        !   so the CFL condition needs to be applied every timestep.
        if (inc_inflows) then
            call dtmax(dtime)
        else
            dtime = dt0
        end if

        ! Adjust timestep to writing and AR intervals
        ! -------------------------------------------
        
        twmaps = .false.
        newars = .false.
        twbfly = .false.
        twavrg = .false.
        twcpoint = .false.
        
        ! Provisional timesteps
        prdtime_maps = dtime + 1.d0
        prdtime_nar = dtime + 1.d0
        prdtime_bfly = dtime + 1.d0
        prdtime_avrg = dtime + 1.d0
        prdtime_cpoint = dtime + 1.d0
        
        if ((wr_b .or. wr_v) .and. tacc_maps + dtime >= dt_write) then
            prdtime_maps = dt_write-tacc_maps
        end if
        if (inc_newar .and. tacc_nar + dtime >= dt_newar) then
            prdtime_nar = dt_newar-tacc_nar
        end if
        if (wr_bfly .and. tacc_bfly + dtime >= dt_bfly) then
            prdtime_bfly = dt_bfly-tacc_bfly
        end if
        if (wr_bfly .and. tacc_avrg + dtime >= dt_avrg) then
            prdtime_avrg = dt_avrg-tacc_avrg
        end if
        if (wr_cp .and. tacc_cpoint + dtime >= dt_cpoint) then
            prdtime_cpoint = dt_cpoint-tacc_cpoint
        end if

        dtime = min(prdtime_maps, prdtime_nar, prdtime_bfly, prdtime_avrg, prdtime_cpoint, dtime)

        if (dtime == prdtime_maps) then
            twmaps = .true.
            tacc_maps = 0.d0
        else
            tacc_maps = tacc_maps + dtime
        end if
            
        if (dtime == prdtime_nar) then
            newars = .true.
            tacc_nar = 0.d0
        else
            tacc_nar = tacc_nar + dtime
        end if
        
        if (dtime == prdtime_bfly) then
            twbfly = .true.
            tacc_bfly = 0.d0
        else
            tacc_bfly = tacc_bfly + dtime
        end if
        
        if (dtime == prdtime_avrg) then
            twavrg = .true.
            tacc_avrg = 0.d0
        else
            tacc_avrg = tacc_avrg + dtime
        end if
        
        if (dtime == prdtime_cpoint) then
            twcpoint = .true.
            tacc_cpoint = 0.d0
        else
            tacc_cpoint = tacc_cpoint + dtime
        end if
        
        if (time+dtime > t_final) then
            dtime = t_final - time
            endloop = .true.
        end if

    end subroutine calc_dt
    
    subroutine dtmax(dtime)
    
        double precision :: maxv, dtime, mphi, mtht, dtimet, dtimep

        dtimet = minval(abs(cfl*dx*rsun/st/vt))
        dtimep = minval(abs(cfl*dphi*rsun*ct/vp))

        dtime = min(dtimet, dtimep)

        ! With an integer number of seconds I avoid extra steps
        ! (writing and placing new ars within a fraction of seconds)
        ! See documentation
        
        dtime = float(floor(dtime))

    end subroutine dtmax
    
end module timestep
