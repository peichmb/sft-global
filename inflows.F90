module inflows
    
    use params
    use grid
    use crank
    use vars
    use derivatives
    use lsflows
    
    implicit none
    
contains

    subroutine add_inflows
    
        integer :: i, j
        double precision, dimension(nphi,nx) :: deriv_aux
        
        ! Compute inflows from magnetic field
        ! -----------------------------------
        b_sm = b
        
        b_sm = abs(b_sm)
        
        if (lfluxinf) call calcf(b_sm)
        
        if (lfiltermode == 1) then
            b_sm = b_sm*st*2.d0/sqrt(3.d0) 
        else if (lfiltermode == 2) then
            b_sm = b_sm*erfilter
        end if
            
        do i = 1, nsteps_sm
            call idiffx(dt_sm,101) ! Term #1 of diffusive smoothing
        !    call idiffp(dt_sm,102) ! Term #2 of diffusive smoothing
        end do

        call sediffp(b_sm,eta,nsteps_sm*dt_sm)
        
        !!!!!!!!!!!!
        ! TEST BLOCK
!         open(736,file = 'b_sm')
!         
!         write(736,*) 0.d0
!         do j = 1, nx
!         write(736,'(<nphi>f30.5)') b_sm(:,j)
!         end do
! 
!         close(736)
         !END OF TEST BLOCK
         !!!!!!!!!!!!!!!!!!

        b_sm = b_sm**inf_b        
        
        ! theta- component of inflow velocity
        call deriv_x_fs(b_sm,deriv_aux)
        vt_inf = -inf_a/rsun*st*deriv_aux

        ! Averaged inflows (phi component vanishes)
        if (latavinf) then
            do j = 1, nx
                vt_inf(:,j) = sum(vt_inf(:,j))/nphi
            end do
        end if

        ! phi- component of inflow velocity
        if (.not. latavinf) then
            call deriv_phi(b_sm,deriv_aux)
            vp_inf = inf_a/rsun/st*deriv_aux
        end if
        
        ! Add inflows to velocity maps              
        if (lfiltermode == 2) then
            vt_inf = vt_inf*erfilter
            vp_inf = vp_inf*erfilter
        end if
        vt = mflow + vt_inf
        if (.not. latavinf) vp = drot + vp_inf

    end subroutine add_inflows
    
    
    subroutine calcf(b_sm)
    
        integer :: i, j
        double precision, dimension(nphi,nx) :: b_sm
        ! Accounts for darkening beyond 200G; see Voegler (2005) Mem. S,A,It. Vol 76, 842
        
        do j = 1, nx
        do i = 1, nphi
        
            if(b_sm(i,j) <= 50.d0) then
                b_sm(i,j) = 2.d-4*b_sm(i,j)
            else if(b_sm(i,j) <= 200.d0) then
                b_sm(i,j) = 0.005d0 + 1.d-4*b_sm(i,j)
            else if(b_sm(i,j) <= 400.d0) then
                b_sm(i,j) = 0.04d0 - 7.5d-5*b_sm(i,j)
            else
                b_sm(i,j) = 0.065d0 - 1.375d-4*b_sm(i,j)
            end if
        
        end do
        end do
        
    end subroutine calcf
        
end module inflows
