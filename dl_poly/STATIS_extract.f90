!       program extract_STATIS.f
!***********************************************************
!       Takes a dl_poly STATIS file and organizes data
!       of interest into columns of an output file.
!***********************************************************
        implicit none

        integer N, i, j, k
        parameter (N=100000)
        character*8 tstep
        character*14 time,P,V,T,E,H,msd
        character*14 Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz
        real d1,d2,d3,d4,d5

        ! Open STATIS file and remove header
        open(1,file='STATIS',status='old',ERR=90)
        read(1,*)
        read(1,*)

        ! Two output files, one containing all thermodynamic
        ! variables (STATIS.dat), the other purely the nine
        ! components of the stress tensor (STRESS.dat)
        ! STATIS.dat: tstep,time,P,V,T,E,H,msd
        ! STRESS.dat: tstep,time,Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz
        ! Energy is usually in eV as specified in the FIELD file
        ! time is in picoseconds, length is in angstroms
        open(2,file='statis.dat')
        open(3,file='stress.dat')

        ! Read in the STATIS file and write to files
        do i=1,N

          ! Read in one block of data
          read(1,*,END=100) tstep,time,d1
          read(1,*) d1,T,E,d4,d5
          read(1,*) d1,d2,d3,d4,H
          read(1,*) d1,d2,d3,d4,d5
          read(1,*) d1,d2,d3,V,d5
          read(1,*) d1,d2,d3,d4,d5
          read(1,*) d1,P,msd,Sxx,Sxy
          read(1,*) Sxz,Syx,Syy,Syz,Szx
          read(1,*) Szy,Szz

          write(2,*) tstep,time,P,V,T,E,H,msd
          write(3,*) tstep,time,Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz
!          namelist/statis/tstep,time,P,V,T,E,H,msd
!          namelist/stress/tstep,time,Sxx,Sxy,Sxz,Syx,Syy,Syz,Szx,Szy,Szz
!          write(2,nml=statis)
!          write(3,nml=stress)

        enddo

 100    continue

        close(1)
        close(2)
        close(3)

 90     continue

        END
