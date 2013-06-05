!       program convert_HISTORY_xyz.f
!***********************************************************
!       Takes a dl_poly HISTORY file and createse a
!       xyz file from it.
!***********************************************************
        implicit none

        integer i, j, k, maxSteps
        parameter (maxSteps=100000)
        integer natom, tstep, levcfg
        real*8 ax, ay, az, bx, by ,bz, cx, cy, cz
        real*8 x, y, z, vx, vy, vz, fx, fy, fz
        character*2 typat, dummy

        open(1,file='HISTORY',status='old',ERR=90)
        open(2,file='HISTORY.xyz')

        ! HISTORY header
        read(1,*)
        read(1,*) levcfg

        ! Read in the HISTORY file and write to xyz file
        do i=1,maxSteps

          read(1,*,END=100) dummy, tstep, natom
          read(1,*) ax, ay, az
          read(1,*) bx, by, bz
          read(1,*) cx, cy ,cz

          write(2,*) natom
          write(2,*) tstep

          do j=1,natom

            read(1,*) typat
            read(1,*) x, y, z
            if (levcfg.gt.0) read(1,*) vx, vy, vz
            if (levcfg.gt.1) read(1,*) fx, fy, fz

 5          format(A2,x,f12.8,x,f12.8,x,f12.8)
            write(2,5) typat, x, y, z

          enddo

        enddo

 100    continue

        close(1)

        close(2)

 90     continue

        END
