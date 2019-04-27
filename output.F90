module output

    use params
    use grid
    use vars
    use iobelda
    
    implicit none
    
    double precision, dimension(:,:), allocatable :: phitav
    character(2) :: suffix
    integer :: avtdim, wrstep = 0
    integer :: fpart = 1, indsave = 1
    
contains

    subroutine init_output()
    
        integer :: i, j
        
        if (verbose) print *,'- Initializing output...'
        
        call int2str(fpart, 2, suffix)
    
        if(wr_bfly) then
            open(21, file=trim(path)//'/bfdiag.'//suffix)
            write(21,*) 'Longitudinal average of the surface field.'
            write(21,*) '-- time [years], B(x) [gauss, x = cos(theta)]'
            write(21,*)
            close(21)
        end if
        
        if(wr_adm) then
            open(22, file=trim(path)//'/axdpm.'//suffix)
            write(22,*) 'Axial dipole moment.'
            write(22,*) '-- time [years], B_p [Gauss]'
            write(22,*)
            close(22)
        end if
        
        !open(211, file=trim(path)//'/totflux.'//suffix)
        !close(211)
        
        avtdim = int(dt_bfly/dt_avrg)
        
        allocate(phitav(avtdim,nx))
        
        !do i = 1, avtdim
        !do j = 1, nx
        !    phitav(i,j) = sum(b(:,j))/float(nphi)
        !end do
        !end do
        
    end subroutine init_output


    subroutine wrmaps(time)

        character(100) :: vp_fname, vt_fname, b_fname
        character(5) :: filenum
        integer :: j
        double precision :: time
        
        call int2str(wrstep,5,filenum)

        if (wr_v) then
        
            if (verbose) print *, 'Writing velocities at t = ',time*secs2years
            
            vp_fname = trim(path)//'/flphi.'//suffix//'.'//filenum
            vt_fname = trim(path)//'/fltht.'//suffix//'.'//filenum
            
            open(11, file = trim(vp_fname))
            open(12, file = trim(vt_fname))
            
            write(11,'(f15.6)') time*secs2years
            write(12,'(f15.6)') time*secs2years
            
            do j = 1, ntht
                write(11,'(<nphi>f20.10)') vp_inf(1:nphi,j)
                write(12,'(<nphi>f20.10)') vt_inf(1:nphi,j)
            end do
            
            close(11)
            close(12)
                            
        end if
        
        if (wr_b) then
        
            if (verbose) print *, 'Writing magnetogram at t = ',time*secs2years, 'years'
        
            b_fname = trim(path)//'/mgram.'//suffix//'.'//filenum
            
            open(13, file = trim(b_fname))
            
            write(13,'(f15.6)') time*secs2years
            
            do j = 1, ntht
                write(13, '(<nphi>f20.10)') b(1:nphi,j)
            end do
            
            close(13)
        
        end if
        
        wrstep = wrstep + 1
        
    end subroutine wrmaps

    
    subroutine wrbfly(time)
    
        double precision :: time
        double precision, dimension(ntht) :: butterfly
        integer :: j
        
        do j = 1, ntht
            butterfly(j) = sum(phitav(:,j))/float(avtdim)
        end do
        
        if (verbose) print *, 'Writing butterfly diagram at t = ', time*secs2years, 'years'
        open(21, position='append', file=trim(path)//'/bfdiag.'//suffix)
        write(21,'(1f20.10,<ntht>f20.10)') time*secs2years, butterfly
        close(21)
        
        !open(211, position='append', file=trim(path)//'/totflux.'//suffix)
        !write(211,*) time*secs2years, sum(b)*(-dx)*dphi*rsuncgs**2/1.e21
        !close(211)
        
    end subroutine wrbfly
    
    
    subroutine wradm(time)
    
        double precision :: time, adm
        
        adm = dx*dphi*sqrt(4.d0/3.d0/pi)*sum(b*ct)

        if (verbose) print *, 'Writing axial dipole moment at t = ', time*secs2years, 'years'
        open(22, position='append', file=trim(path)//'/axdpm.'//suffix)
        write(22,'(2f20.10)') time*secs2years, adm
        close(22)                
    
    end subroutine wradm

    
    subroutine saveforav()
    
        integer :: j
    
        do j = 1, nx
            phitav(indsave,j) = sum(b(:,j))/float(nphi)
        end do
        
        indsave = indsave+1
        
        if(indsave > avtdim) indsave = 1
            
    end subroutine saveforav
    
end module output
