!       program supercell.f90
!***********************************************************
!       Create a 2x2x2 supercell xyz file from an
!       original wrapped xyz file.
!***********************************************************
        implicit none
        
        integer maxSteps, i, j, k
        parameter (maxSteps=1000000)
        integer natom, tstep
        real x, y, z, ax, ay, az
        character*8 typat
        character*128 fin

        ! Retrieve required user input
        write(*,*) 'Name of wrapped xyz file'
        read(*,*) fin
        write(*,*) 'Lattice constants (a,b,c) (angstroms):'
        read(*,*) ax, ay, az

        ! Open input and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='SUPERCELL.xyz')

        do i=1,maxSteps

          ! Read and write header information
          read(1,*,END=100) natom
          read(1,*) tstep
          write(2,*) natom*2*2*2
          write(2,*) tstep

          do j=1,natom

            ! Read in one coordinate from original xyz file
            read(1,*) typat, x, y, z

            ! Write supercell coordinates
            write(2,*) typat, x,    y,    z
            write(2,*) typat, x+ax, y,    z
            write(2,*) typat, x+ax, y+ay, z
            write(2,*) typat, x+ax, y,    z+az
            write(2,*) typat, x+ax, y+ay, z+az
            write(2,*) typat, x,    y+ay, z
            write(2,*) typat, x,    y+ay, z+az
            write(2,*) typat, x,    y,    z+az

          enddo

        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
