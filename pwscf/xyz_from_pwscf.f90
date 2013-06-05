!       program xyz_from_pwscf.f
!***********************************************************
!       Create an xyz file from pwscf output files
!       must have a temp.xyz file present (made with
!       xyz_from_pwscf.py
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, i, j, k
        parameter (n=10000000,m=100000)
        integer natom
        real x, y, z, a, b, c
        character*8 typat

        ! Open the input and output files
        open(1,file='temp.xyz',status='old',ERR=90)
        open(2,file='TRAJEC.xyz')

        ! Retrieve necessary input
        write(*,*) 'Number of atoms'
        read(*,*) natom
        write(*,*) 'Lattice constants (in angstroms); a,b,c'
        read(*,*) a, b, c

        ! Loop over all timesteps
        do i=1,n

          ! Write natom and timestep to xyz file
          write(2,*) natom
          write(2,*) i

          ! Loop over natom
          do j=1,natom

            ! Read in a coord from temp.xyz
            read(1,*,END=100) typat, x, y, z
            
            ! Convert coords to angstroms
            x = x*a
            y = y*b
            z = z*c

            ! Write converted coords to TRAJEC.xyz
            write(2,*) typat, x, y, z

          enddo

        enddo        

 100    continue

        close(1)

        close(2)

 90     continue

        END
