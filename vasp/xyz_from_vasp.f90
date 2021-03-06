!       program xyz_from_vasp.f
!***********************************************************
!       Create an xyz file from vasp output files
!       Made from XDATCAR file
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer maxSteps, maxAtoms
        parameter (maxSteps=10000000,maxAtoms=100000)
        integer i, j
        integer natom
        real x, y, z, alat, a, b, c, d0, d1, d2
        character*2 typat

        ! Open the input and output files
        open(1,file='XDATCAR',status='old',ERR=90)
        open(2,file='POSCAR',status='old',ERR=80)
        open(3,file='TRAJEC.xyz')

        ! Get cell and natom from POSCAR file
        read(2,*)
        read(2,*) alat
        read(2,*) a, d1, d2
        read(2,*) d0, b, d2
        read(2,*) d0, d1, c
        read(2,*) natom
        close(2)

        a = a*alat
        b = b*alat
        c = c*alat

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

            ! Read in
            read(1,*,END=100,ERR=111) x, y, z
            
            ! Convert coords from reduced to angstroms
            x = x*a
            y = y*b
            z = z*c

            ! Write converted coords to TRAJEC.xyz
            write(3,*) typat, x, y, z

          enddo

        enddo

 111    continue

        write(*,*) "ERROR ON LINE: ", i, j

 100    continue

        close(3)

 80     continue

        close(1)

 90     continue

        END
