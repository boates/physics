!       program xyz_from_vasp_general.f90
!***********************************************************
!       Create an xyz file from vasp output files
!       Made from XDATCAR file
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer maxSteps, maxAtoms, i, j, natom
        parameter (maxSteps=10000000,maxAtoms=100000)
        real x, y, z, x_r, y_r, z_r
        real alat, ax, ay, az, bx, by, bz, cx, cy, cz
        character*2 typat

        ! Open the input and output files
        open(1,file='XDATCAR',status='old',ERR=90)
        open(2,file='POSCAR',status='old',ERR=80)
        open(3,file='TRAJEC.xyz')

        ! Get cell and natom from POSCAR file
        read(2,*)
        read(2,*) alat
        read(2,*) ax, ay, az
        read(2,*) bx, by, bz
        read(2,*) cx, cy, cz
        read(2,*) natom
        close(2)

        ax = ax*alat
        ay = ay*alat
        az = az*alat
        bx = bx*alat
        by = by*alat
        bz = bz*alat
        cx = cx*alat
        cy = cy*alat
        cz = cz*alat

        ! Read off XDATCAR header and get typat
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) typat

        typat = 'N'

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
            read(1,*,END=100) x_r, y_r, z_r
            
            ! Convert coords from reduced to angstroms
            x = x_r*ax + y_r*bx + z_r*cx
            y = x_r*ay + y_r*by + z_r*cy
            z = x_r*az + y_r*bz + z_r*cz

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
