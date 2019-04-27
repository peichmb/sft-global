module initialize

    use params
    use random
    use grid
    use sources
    use derivatives
    use vars
    use initconds
    use lsflows
    use rk4
    use timestep
    use crank
    use inflows
    use update
    use output
    use cpoint

    implicit none
    
contains

    subroutine init(time, step, res)
    
        integer :: step, i, rt
        double precision :: time
        logical :: res
    
        if(verbose) then
            print *
            print *, 'Initializing modules and setting initial conditions.'
        end if

        if(res) then
            ! In this case params, vars and grid are initialized in the module cpoint
            call resume(time, step)
            call init_random()
            call init_lsflows()
        else
            call init_params(res)
            call init_random()
            call init_grid()
            call init_vars()
            call init_lsflows()
            call set_initconds()
            call init_output()
        end if
        
        if (verbose) call show_params()

        !print *, inc_newar 
        !print *, dt_newar
        if (inc_newar) call init_sources()
        !print *, dt_newar
        
        ! Set velocities
        vp = drot
        vt = mflow
        
        call init_timestep()
        call init_derivatives()        
        call init_rk4()
        call init_crank()

        if (inc_inflows) call add_inflows()

    end subroutine init
    
end module initialize
