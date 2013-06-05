!       program PFB_from_vasp.f90
!***********************************************************
!       Make position_force_bonding file from cbn file and
!       FORCES.xyz file. FORCES.xyz is created by
!       /home/boates/software/vasp/FORCES_from_OUTCAR.py
!       and is in eV/angst. PFB file will be written in
!       angst and eV/angst.
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, nn
        real x, y, z, fx, fy, fz
        character*8 typat, d1, d2, d3
        character*128 fcbn

        ! Retrieve user input
        write(*,*) 'Name of cbn file'
        read(*,*) fcbn
        write(*,*) 'Assuming FORCES.xyz file is present and of same length as:',fcbn
        
        ! Open input and output files
        open(1,file=fcbn,status='old',ERR=90)
        open(2,file='FORCES.xyz',status='old',ERR=90)
        open(3,file='PFB.dat')

        ! Read cbn file header
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

 5      format(f10.6,x,f10.6,x,f10.6,x,f10.6,x,f10.6,x,f10.6,x,i3)

        ! Begin routine
        do i=1,maxSteps

          do j=1,natom

            ! Read coords and forces
            read(1,*,END=100) typat, x, y, z, nn
            read(2,*) fx, fy, fz

            ! convert coords to angstroms and write to file
            x = x * 0.529177
            y = y * 0.529177
            z = z * 0.529177
            write(3,5) x, y, z, fx, fy, fz, nn

          enddo
          
        enddo

 100    continue

        close(1)
        close(2)
        close(3)

 90     continue

        END
