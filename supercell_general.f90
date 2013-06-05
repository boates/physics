!       program supercell_general.f90
!***********************************************************
!       Create a 2x2x2 supercell xyz file from an
!       original wrapped xyz file.
!***********************************************************
        implicit none
        
        integer maxSteps, i, j, k, l, m
        parameter (maxSteps=1000000)
        integer natom, tstep, nx, ny, nz
        real x, y, z, alat, ax, ay, az, bx, by, bz, cx, cy, cz
        character*8 typat
        character*128 fin, car

        ! Retrieve required user input
        write(*,*) 'Name of wrapped xyz file'
        read(*,*) fin
        write(*,*) 'Name of POSCAR or CONTCAR file'
        read(*,*) car
        write(*,*) 'Replications in x,y,z:'
        read(*,*) nx, ny, nz

        ! Open input and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file=car,status='old',ERR=90)
        open(3,file='SUPERCELL.xyz')

        ! Read lattice vectors from CONTCAR file
        read(2,*)
        read(2,*) alat
        read(2,*) ax, ay, az
        read(2,*) bx, by, bz
        read(2,*) cx, cy, cz
        close(2)
        
        ax = alat*ax
        ay = alat*ay
        az = alat*az
        bx = alat*bx
        by = alat*by
        bz = alat*bz
        cx = alat*cx
        cy = alat*cy
        cz = alat*cz

        write(*,*) 'New lattice vectors (angstroms)'
        write(*,*) 'a = ',nx*ax,nx*ay,nx*az
        write(*,*) 'b = ',ny*bx,ny*by,ny*bz
        write(*,*) 'c = ',nz*cx,nz*cy,nz*cz

        do i=1,maxSteps

          ! Read and write header information
          read(1,*,END=100) natom
          read(1,*) tstep
          write(3,*) natom*nx*ny*nz ! natom*(1 + 7*(nx-1)*(ny-1)*(nz-1))
          write(3,*) tstep

          do j=1,natom

            ! Read in one coordinate from original xyz file
            read(1,*)  typat, x, y, z
            write(3,*) typat, x, y, z

            ! Write supercell coordinates
            do k=1,nx-1
              write(3,*) typat, x+k*ax, y+k*ay, z+k*az
            enddo

            do l=1,ny-1
              write(3,*) typat, x+l*bx, y+l*by, z+l*bz
            enddo

            do m=1,nz-1
              write(3,*) typat, x+m*cx, y+m*cy, z+m*cz
            enddo

            do k=1,nx-1
              do l=1,ny-1
                write(3,*) typat, x+k*ax+l*bx, y+k*ay+l*by, z+k*az+l*bz
              enddo
            enddo

            do k=1,nx-1
              do m=1,nz-1
                write(3,*) typat, x+k*ax+m*cx, y+k*ay+m*cy, z+k*az+m*cz
              enddo
            enddo

            do l=1,ny-1
              do m=1,nz-1
                write(3,*) typat, x+l*bx+m*cx, y+l*by+m*cy, z+l*bz+m*cz
              enddo
            enddo

            do k=1,nx-1
              do l=1,ny-1
                do m=1,nz-1
                  write(3,*) typat, x+k*ax+l*bx+m*cx, y+k*ay+l*by+m*cy, z+k*az+l*bz+m*cz
                enddo
              enddo
            enddo

          enddo

        enddo

 100    continue

        close(3)

 90     continue

        close(1)

        END
