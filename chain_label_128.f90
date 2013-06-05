!       program chain_label.f90
!***********************************************************
!       Read in a CHAIN.cbn file created by chain_angles.py
!       and write an dat file for VMD giving the color 
!       indicator for each atom at each timestep.
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, tstep, nn(maxAtoms), label(maxAtoms)
        real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
        character*2 typat
        character*128 fin,d1,d2,d3,d4,d5

        ! Retrieve user input
        write(*,*) 'Name of cbn file'
        read(*,*) fin

        ! Open input and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='CHAIN_LABEL.dat')

        ! Read cbn file header
        read(1,*,END=100) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) d1,d2,d3,natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        if (natom.ne.128) then
          write(*,*) 'WARNING - chain_label_128.x assumes 128 atom for formatting'
          write(*,*) 'FOUND: natom =',natom,' exiting...'
          goto 100
        endif

        do i=1,maxSteps

          do j=1,natom

            read(1,*,END=100) typat, X(j), Y(j), Z(j), nn(j)

            ! Zero for molecule, One for chain
            if (nn(j).ge.0) then
              label(j) = 0
            endif
            if (nn(j).eq.-1) then
              label(j) = 1
            endif

          enddo
          
 5        format(128i3)
          write(2,5) (label(k),k=1,natom)

        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
