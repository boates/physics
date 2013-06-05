!       program two_phase_rename_half_xyz.f90
!****************************************************************
!       Rename the first 64 atoms in a 128 atom xyz file to
!       something else. Purpose is to visualize two phase
!        simulations
!****************************************************************
        implicit none
        
        integer maxSteps, maxAtoms, maxBins, i, j, k
        parameter (maxSteps=100000,maxAtoms=128, maxBins=1000)
        integer natom, tstep
        real*4 x, y, z
        character*2 typat
        character*128 fxyz

        ! Retrieve user input
        write(*,*) 'Name of .xyz file'
        read(*,*) fxyz

        ! Open .xyz and output files
        open(1,file=fxyz,status='old',ERR=90)
        open(2,file='TWO_PHASE.xyz')

        ! Perform algorithm
        do i=1,maxSteps

          ! Read & write natom and tstep
          read(1,*,END=100) natom
          read(1,*) tstep
          write(2,*) natom
          write(2,*) tstep

          ! Read and write coordinates
          do j=1,natom/2
            read(1,*) typat, x, y, z
            write(2,*) 'C', x, y, z
          enddo

          do j=1,natom - natom/2
            read(1,*) typat, x, y, z
            write(2,*) typat, x, y, z
          enddo

        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
