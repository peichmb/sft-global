program sft

    use params
    use initialize
    use sources
    use timestep
    use inflows
    use update
    use output
    use cpoint
    
    implicit none
    
    integer :: step = 0
    double precision :: time = 0.d0, dtime
    ! Time to write maps, place new ARs or add an entry to the bfly. diagram,
    !         average, write checkpoint 
    logical :: twmaps = .false., newars = .false., twbfly = .false., &
             & twavrg = .false., twcpoint = .false.
    ! Time to end the simulation
    logical :: endloop = .false., res = .false.
    
    print *
    print *, '-----------------------------------'
    print *, '--- SURFACE FLUX TRANSPORT CODE ---'
    print *, '-----------------------------------'
    
    ! Are we resuming a simulation?
    if ( IARGC() > 0 ) then
        res = .true.
    end if
    
    ! Initialize modules, set initial conditions and calculate initial timestep
    call init(time, step, res)
    
    ! Write initial state
    !if (.not. res) then
        call wrmaps(time)
!        call wradm(time)
        
    !end if
    
    ! Main loop
    print *, 'Starting simulation at t = ', time*secs2years,' years'
    do while(.not.endloop)
        
        ! Calculate timestep
        call calc_dt(dtime, twmaps, newars, twbfly, twavrg, twcpoint, &
                   & time, t_final, endloop)
        if (verbose) write(*,'(a9,i10,a18,f7.3,a12,f9.1,a2)') &
         & '| Step = ', step, ' | time [years] = ', time*secs2years, ' | dt [s] = ', dtime,' |'
        
        ! Update solution to time+dtime
        call upsol(dtime)
        time = time + dtime
        step = step + 1
        if (step == maxsteps) endloop = .true.
        
        ! Place new ARs if it is time
        if (newars) call placenewars(time)

        ! Compute inflows from magnetic field
        if (inc_inflows) call add_inflows()

        ! Write solution if it is time
        if (twmaps) call wrmaps(time)
        !call wrmaps(time)
        ! Average
        if (twavrg) call saveforav()
        
        ! Add entry to butterfly diagram if it is time
        if (twbfly) then
            if (wr_bfly) call wrbfly(time)
            if (wr_adm) call wradm(time)
        end if
                   
        ! Write check point if it is time
        if (twcpoint) call checkpoint(time,step)                   
        
    end do
    
    print *, 'Simulation ended at time = ',time*secs2years,' years.'
    print *, 'Number of steps: ', step
    print *
    
    if (armode==2) close(321)
    
end program sft
