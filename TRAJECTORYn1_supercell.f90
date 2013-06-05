!       program TRAJECTORYn1_supercell.f90
!***********************************************************
!       Create a nx X ny X nz supercell xyz file from an
!       original TRAJECTORYn1 file
!***********************************************************
        implicit none
        
        integer maxSteps, i, j, k, l, m
        parameter (maxSteps=1000000)
        integer natom, tstep, nx, ny, nz
        real*8 x, y, z, alat, ax, ay, az, bx, by, bz, cx, cy, cz
        real*8 fx, fy, fz
        character*128 fin, car

        ! Retrieve required user input
        write(*,*) 'Name of TRAJECTORYn1 file'
        read(*,*) fin
        write(*,*) 'Total number of atoms'
        read(*,*) natom
        write(*,*) 'Name of POSCAR or CONTCAR file'
        read(*,*) car
        write(*,*) 'Replications in x,y,z:'
        read(*,*) nx, ny, nz

        ! Open input and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file=car,status='old',ERR=90)
        open(3,file='TRAJECTORYn1.supercell')

        ! Read lattice vectors from POSCAR/CONTCAR file
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

 4      format(i2)
 5      format(f12.8,x,f12.8,x,f12.8,x,f12.8,x,f12.8,x,f12.8)

!        write(*,*) 'New lattice vectors (angstroms)'
!        write(*,*) 'a = ',nx*ax,nx*ay,nx*az
!        write(*,*) 'b = ',ny*bx,ny*by,ny*bz
!        write(*,*) 'c = ',nz*cx,nz*cy,nz*cz

        do i=1,maxSteps

          ! Read and write header information
!          read(1,*,END=100) natom
          read(1,*,END=100) tstep
!          write(3,*) natom*nx*ny*nz ! natom*(1 + 7*(nx-1)*(ny-1)*(nz-1))
          write(3,4) tstep

          do j=1,natom

            ! Read in one coordinate from original xyz file
            read(1,*,END=100) x, y, z, fx, fy, fz
            write(3,5) x, y, z, fx, fy, fz

            ! Write supercell coordinates
            do k=1,nx-1
              write(3,5) x+k*ax, y+k*ay, z+k*az, fx, fy, fz
            enddo

            do l=1,ny-1
              write(3,5) x+l*bx, y+l*by, z+l*bz, fx, fy, fz
            enddo

            do m=1,nz-1
              write(3,5) x+m*cx, y+m*cy, z+m*cz, fx, fy, fz
            enddo

            do k=1,nx-1
              do l=1,ny-1
                write(3,5) x+k*ax+l*bx, y+k*ay+l*by, z+k*az+l*bz, fx, fy, fz
              enddo
            enddo

            do k=1,nx-1
              do m=1,nz-1
                write(3,5) x+k*ax+m*cx, y+k*ay+m*cy, z+k*az+m*cz, fx, fy, fz
              enddo
            enddo

            do l=1,ny-1
              do m=1,nz-1
                write(3,5) x+l*bx+m*cx, y+l*by+m*cy, z+l*bz+m*cz, fx, fy, fz
              enddo
            enddo

            do k=1,nx-1
              do l=1,ny-1
                do m=1,nz-1
                  write(3,5) x+k*ax+l*bx+m*cx, y+k*ay+l*by+m*cy, z+k*az+l*bz+m*cz, fx, fy, fz
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
