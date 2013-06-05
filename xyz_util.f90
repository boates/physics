!       program xyz_util.f
!***********************************************************
!       Perform operations on an xyz file (TRAJEC.xyz)
!       1) replace negative coords with PBC displacements
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, i, j, k
        parameter (n=1000000,m=10000)
        integer natom, timestep
        real*4 x, y, z
        real*8 ax, ay, az
        character*2 typat

        ! Get the lattice constants from the user
        write(*,*) 'Lattice constants: ax, ay, az'
        read(*,*) ax, ay, az

        ! xyz file assumed to be named TRAJEC.xyz
        ! new file will be named new.TRAJEC.xyz
        open(1,file='TRAJEC.xyz',status='old',ERR=90)
        open(2,file='new.TRAJEC.xyz')

        ! Read in the value of natom and timestep
        read(1,*) natom
        read(1,*) timestep

        ! Loop over the rest of the file
        do i=1,n

          ! Write natom and timestep to new xyz file
          write(2,*) natom
          write(2,*) i

          do j=1,natom

            ! Read in coordinates from xyz file
            read(1,*,END=100) typat, x, y, z

            ! Checks for negative coordinates
            if (x.lt.0.0) then
              x = x + ax
            endif
            if (y.lt.0.0) then
              y = y + ay
            endif
            if (z.lt.0.0) then
              z = z + az
            endif

            write(2,*) typat, x, y, z

          enddo

        ! Read off the new natom and t
        read(1,*,END=100) natom
        read(1,*) timestep

        enddo

 100    continue

        close(1)

 90     continue

        END
