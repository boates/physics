!       program xyz_from_cbn.f90
!***********************************************************
!       Convert cbn file back to xyz file.
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom
        real x, y, z
        character*8 typat,d1,d2,d3
        character*128 fin

        ! Retrieve user input
        write(6,*) 'Name of .cnn or .cbn file'
        read(5,*) fin

        ! Open input and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='TRAJEC.xyz')

        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) d1, d2, d3, natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        do i=1,maxSteps

          read(1,*,END=100) typat, x, y, z
          write(2,*) natom
          write(2,*) i
          write(2,*) typat, x*0.529177, y*0.529177, z*0.529177

          do j=1,natom-1

            read(1,*,END=100) typat, x, y, z
            write(2,*) typat, x*0.529177, y*0.529177, z*0.529177

          enddo
          
        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
