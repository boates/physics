!       program PBC_extend_x.f90
!***********************************************************
!       Create a periodically repeated xyz file in +/- x
!***********************************************************
        implicit none
        
        integer maxSteps, i, j, k
        parameter (maxSteps=1000000)
        integer natom, tstep
        real x, y, z, ax
        character*8 typat
        character*128 fin

        ! Retrieve required user input
        write(*,*) 'Name of wrapped xyz file'
        read(*,*) fin
        write(*,*) 'Lattice constant (ax) (angstroms):'
        read(*,*) ax

        ! Open input and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='PBC_EXTEND_X.xyz')

        do i=1,maxSteps

          ! Read and write header information
          read(1,*,END=100) natom
          read(1,*) tstep
          write(2,*) natom + natom*2
          write(2,*) tstep

          do j=1,natom

            ! Read in one coordinate from original xyz file
            read(1,*) typat, x, y, z

            ! Write extended coordinates
            write(2,*) typat, x,    y,    z
            write(2,*) typat, x+ax, y,    z
            write(2,*) typat, x-ax, y,    z

          enddo

        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
