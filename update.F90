module update
    
    use params    
    use rk4
    use crank
    use vars
    use grid
    use derivatives
    use timestep
    
    implicit none
    
contains

    subroutine upsol(dtime)
        
        double precision :: dtime
        
        call advection(dtime)
        call diffusion(dtime)
    
    end subroutine upsol

    subroutine advection(dtime)
    
        double precision :: dtime
    
        call rkstep(dtime,b,1) ! Term #1 of SFTE (x- advection)
        call rkstep(dtime,b,2) ! Term #2 of SFTE (phi - advection)
    
    end subroutine advection
    
    subroutine diffusion(dtime)
    
        double precision :: dtime

        call idiffx(dtime,3) ! Term #3 of SFTE (x- diffusion)
        !call idiffp(dtime,4) ! Term #4 of SFTE (phi- diffusion)
        call sediffp(b, eta, dtime)

        !!!!!!!!!! WRONG! dt can also change bc of the writing parameters !! if (inc_inflows) dt0 = dtime ! Necessary to recalculate alpha and beta in the next time-step
        dt_prev = dtime
        
    end subroutine diffusion   
    
end module update
