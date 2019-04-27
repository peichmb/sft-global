module cpoint

    use random
    use params
    use grid
    use iobelda
    use vars
    use output
    use timestep
    
    implicit none
    
    integer :: cpcount = 0
    
contains

    subroutine checkpoint(time, step)
    
        integer :: step, j
        double precision :: time
        character(100) :: fname
        character(5) :: filenum

        cpcount = cpcount + 1
                
        call int2str(cpcount,5,filenum)
        
        fname=trim(path)//'/checkpoint.'//filenum
        
        open(1010, file = trim(fname) )
        
        write(1010, '(f30.5,7i15)') time, step, rnd_rcount, rnd_ncount, cpcount, wrstep, fpart, indsave
        write(1010, '(5f30.5)') tacc_nar, tacc_bfly, tacc_maps, tacc_avrg, tacc_cpoint        
        write(1010, '(2i15,4f30.5,2i15)') nphi, ntht, cfl, pcbnd, etakm2s, t_final*secs2years, maxsteps, rseed
        write(1010, '(2l5,f30.5,i15,2f30.5,l5)') inc_inflows, latavinf, fwhm_sm, nsteps_sm, inf_a, inf_b, lfluxinf
        write(1010, '(i15,2f30.5)') lfiltermode, latfilter, wfilter
        write(1010, '(a100)') scenario
        write(1010, '(3f30.5)') icondp1, icondp2, icondp3
        write(1010, '(l5,i15)') inc_newar, armode
        write(1010, '(a100)') inputrecord
        write(1010, '(3f30.5,l5,8f30.5)') cclength*secs2years, ccoverlap*secs2years, tau0*secs2years, normstpol, joysfac, angscat, actfac, alatu, alatscatu, alatl, alatscatl, bmax
        write(1010, '(l5)') verbose
        write(1010, '(a100)') path
        write(1010, '(2l5,f30.5,2l5,f30.5,l5,f30.5)') wr_v, wr_b, dt_write*secs2days, wr_bfly, wr_adm, dt_bfly*secs2days, wr_cp, dt_cpoint*secs2years
        
        do j = 1, ntht
            write(1010, '(<nphi>f)') b(1:nphi,j)
        end do
        do j = 1, ntht
            write(1010, '(<avtdim>f)') phitav(1:avtdim,j)
        end do
        
        close(1010)
    
        if (verbose) print *, 'Check point at t = ', time*secs2years, 'years'
    
    end subroutine checkpoint
    
    subroutine resume(time, step)
    
        integer :: rnd_totcount, step, fn, j
        double precision :: time
        character(100) :: path_cp
        character(100) :: fname
        character(5) :: filenum

        if (IARGC() < 2) then
            print *, '-- checkpoint.F90; resume --'
            print *, 'Missing arguments.'
            stop
        end if

        call GETARG(1, path_cp)
        call GETARG(2, filenum)
        read(filenum,*) fn
        call int2str(fn,5,filenum)
        fname = trim(path_cp)//'/checkpoint.'//filenum
        
        open(1010, file = trim(fname) )
        
        read(1010, '(f30.5,7i15)') time, step, rnd_rcount, rnd_ncount, cpcount, wrstep, fpart, indsave
                !print *, time, step, rnd_rcount, rnd_ncount, cpcount, wrstep, fpart
        read(1010, '(5f30.5)') tacc_nar, tacc_bfly, tacc_maps, tacc_avrg, tacc_cpoint        
                !print *, tacc_nar, tacc_bfly, tacc_maps, tacc_avrg, tacc_cpoint        
        read(1010, '(2i15,4f30.5,2i15)') nphi, ntht, cfl, pcbnd, etakm2s, t_final, maxsteps, rseed
                !print *, nphi, ntht, cfl, pcbnd, etakm2s, t_final, maxsteps, rseed
        read(1010, '(2l5,f30.5,i15,2f30.5,l5)') inc_inflows, latavinf, fwhm_sm, nsteps_sm, inf_a, inf_b, lfluxinf
                !print *, inc_inflows, latavinf, fwhm_sm, eta_sm, nsteps_sm, inf_a, inf_b, b_thr
        read(1010, '(i15,2f30.5)') lfiltermode, latfilter, wfilter
                !print *, lfiltermode, latfilter, wfilter
        read(1010, '(a100)') scenario
                !print *, scenario
        read(1010, '(3f30.5)') icondp1, icondp2, icondp3
                !print *, icondp1, icondp2, icondp3
        read(1010, '(l5,i15)') inc_newar, armode
                !print *, inc_newar, armode
        read(1010, '(a100)') inputrecord
                !print *, inputrecord
        read(1010, '(3f30.5,l5,8f30.5)') cclength, ccoverlap, tau0, normstpol, joysfac, angscat, actfac, alatu, alatscatu, alatl, alatscatl, bmax
                !print *, cclength, ccoverlap, tau0, normstpol, joysfac, angscat, actfac, alatu, alatscatu, alatl, alatscatl, bmax
        read(1010, '(l5)') verbose
                !print *, verbose
        read(1010, '(a100)') path
                !print *, path
        read(1010, '(2l5,f30.5,2l5,f30.5,l5,f30.5)') wr_v, wr_b, dt_write, wr_bfly, wr_adm, dt_bfly, wr_cp, dt_cpoint
                !print *, wr_v, wr_b, dt_write, wr_bfly, wr_adm, dt_bfly, wr_cp, dt_cpoint

        if(IARGC()==3) call read_outparams(time)

        call init_params(.true.)

        call init_grid()
        call init_vars()
        
        fpart = fpart + 1
        
        call init_output()
        
        print *, nphi
        print *, dt_cpoint
        print *, shape(b) 
        
        do j = 1, ntht
            read(1010, '(<nphi>f)') b(1:nphi,j)
        end do
        do j = 1, ntht
            read(1010, '(<avtdim>f)') phitav(1:avtdim,j)
        end do
        
        close(1010)

            
    end subroutine resume
    
    subroutine read_outparams(time)
    
        integer :: line
        double precision :: t0, time
        character(5) :: t0str
        
        !tacc_nar = 0.d0
        tacc_bfly = 0.d0
        tacc_maps = 0.d0
        tacc_avrg = 0.d0
        tacc_cpoint = 0.d0
        
        open(10,file='parameters')
        do line = 1,89
            read(10,*)
        end do
        
        call GETARG(3,t0str)
        read(t0str,*) t0

        t_final = time/3600.d0/24.d0/365.d0 + t0
        
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
        
    end subroutine read_outparams
    
end module cpoint

