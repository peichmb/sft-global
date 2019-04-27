module random
    
    use params
    
    implicit none
    
    include "interfaces.F90"
    
    integer :: rnd_rcount = 0, rnd_ncount = 0
    double precision :: dummy = 1.d0
    
contains

    subroutine init_random()
    
        integer :: i
    
        if (verbose) print *,'- Initializing random...'
    
        call ZBQLINI(rseed)

        do i = 1, rnd_rcount
            dummy = ZBQLU01(dummy)
        end do
        
        do i = 1, rnd_ncount
            dummy = ZBQLNOR(0.d0,0.d0)
        end do
        
    end subroutine init_random

    function rnd_u01()

        double precision :: rnd_u01
        rnd_rcount = rnd_rcount+1
        rnd_u01 = ZBQLU01(dummy)
        
    end function rnd_u01
    
    function rnd_nor(mu, sigma)
    
        double precision :: mu, sigma, rnd_nor
        
        rnd_ncount = rnd_ncount + 1
        rnd_nor = ZBQLNOR(mu, sigma)
        
    end function rnd_nor
    
    function rndtot()
    
        integer :: rndtot
        
        rndtot = rnd_rcount + rnd_ncount/2*2+2
    
    end function rndtot
    
end module random

