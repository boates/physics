!       program xred_from_vasp.f90
!***********************************************************
!       Create an xyz file in reduced coords from XDATCAR
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer maxSteps, maxAtoms, i, j, natom
        parameter (maxSteps=10000000,maxAtoms=100000)
        real x, y, z, x_r, y_r, z_r
        character*2 typat

        ! Open the input and output files
        open(1,file='XDATCAR',status='old',ERR=90)
        open(3,file='XRED.xyz')

        ! Read off XDATCAR header and get typat
        read(1,*) natom
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) typat

        ! Loop over timesteps
        do i=1,maxSteps

          ! Read off the whitespace line
          read(1,*,END=100)

          ! Write natom and timestep to xyz file
          write(3,*) natom
          write(3,*) i

          ! Loop over natom
          do j=1,natom

            ! Read in reduced coordinates
            read(1,*,END=100) x, y, z
            
            ! Write converted coords to TRAJEC.xyz
            write(3,*) typat, x, y, z

          enddo

        enddo

 100    continue

        close(3)

 80     continue

        close(1)

 90     continue

        END
