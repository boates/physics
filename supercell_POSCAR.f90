!       program supercell_POSCAR.f90
!***********************************************************
!       Create a nx x ny x nz supercell POSCAR file from an
!       original wrapped POSCAR/CONTCAR file.
!***********************************************************
        implicit none
        
        integer maxSteps, i, j, k, l, m
        parameter (maxSteps=1000000)
        integer natom, tstep, nx, ny, nz
        real x, y, z, alat, ax, ay, az, bx, by, bz, cx, cy, cz
        character*8 typat
        character*128 car, header

        ! Retrieve required user input
        write(*,*) 'Name of POSCAR or CONTCAR file'
        read(*,*) car
        write(*,*) 'Replications in x,y,z:'
        read(*,*) nx, ny, nz

        ! Open input and output files
        open(1,file=car,status='old',ERR=90)
        open(2,file='supercell.POSCAR')

        ! Read and write header information
        read(1,*) header
        read(1,*) alat
        read(1,*) ax, ay, az
        read(1,*) bx, by, bz
        read(1,*) cx, cy, cz
        read(1,*) natom
        read(1,*)

 5      format(xxxx,f12.8,x,f12.8,x,f12.8)
        write(2,*) 'POSCAR created by supercell_POSCAR.x'
        write(2,*) '    1.0000000000'
        write(2,5) nx*ax*alat, nx*ay*alat, nx*az*alat
        write(2,5) ny*bx*alat, ny*by*alat, ny*bz*alat
        write(2,5) nz*cx*alat, nz*cy*alat, nz*cz*alat
        write(2,'(x,i6)') natom*nx*ny*nz
        write(2,*) 'Direct'

 6      format(x,f12.8,x,f12.8,x,f12.8)
        do j=1,natom

          ! Read in one coordinate from original POSCAR/CONTCAR file
          read(1,*)  x, y, z
          write(2,6) x/nx, y/ny, z/nz

          ! Write supercell coordinates
          do k=1,nx-1
            write(2,6) (x+k)/nx, y/ny, z/nz
          enddo

          do l=1,ny-1
            write(2,6) x/nx, (y+l)/ny, z/nz
          enddo

          do m=1,nz-1
            write(2,6) x/nx, y/ny, (z+m)/nz
          enddo

          do k=1,nx-1
            do l=1,ny-1
              write(2,6) (x+k)/nx, (y+l)/ny, z/nz
            enddo
          enddo

          do k=1,nx-1
            do m=1,nz-1
              write(2,6) (x+k)/nx, y/ny, (z+m)/nz
            enddo
          enddo

          do l=1,ny-1
            do m=1,nz-1
              write(2,6) x/nx, (y+l)/ny, (z+m)/nz
            enddo
          enddo

          do k=1,nx-1
            do l=1,ny-1
              do m=1,nz-1
                write(2,6) (x+k)/nx, (y+l)/ny, (z+m)/nz
              enddo
            enddo
          enddo

        enddo

 100    continue

        close(2)

 90     continue

        close(1)

        END
