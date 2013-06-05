!       program unwrap_PBC.f
!******************************************************************
!       Unwrap periodically wrapped coordinates from an xyz file
!       Output the unwrapped coords to a new xyz file
!
!       NOTE: This algorithm assumes a reasonably small stepsize
!             with respect to the lattice constants. Is partially
!             taken from F.26 from Allen-Tildesley.
!******************************************************************
        implicit none

        ! Declare all required variables
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, i, j, k
        parameter (n=100000,m=512)
        integer natoms(n), timesteps(n)
        character(2) typat(n,m)
        character(128) fin
        real(4) Rx(n,m),Ry(n,m),Rz(n,m)
        real(4) dx(n,m),dy(n,m),dz(n,m),dRx,dRy,dRz
        real(4) a, ax, ay, az

        ! Ask user for required input, which can be given via
        ! an input file (still echoes the requests though...)
        write(*,*) "Name of xyz file with wrapped coordinates"
        read(*,*) fin
        write(*,*) "Lattice constants (ax,ay,az)"
        read(*,*) ax, ay, az

        ! Assume cubic cell:
!        ax = a
!        ay = a
!        az = a

        ! Open the chosen input file in "unit" 1
        ! Which must be simply 3 columns of coordinates
        open(1,file=fin,status='old',ERR=90)
        ! Open an output file of the chosen name in "unit" 2
        open(2,file='unwrapped.xyz')

        ! Read in the wrapped xyz file data
        ! Convert to reduced coordinates and shift origin
        do i=1,n
          read(1,*,END=100) natoms(i)
          read(1,*) timesteps(i)
          do j=1,natoms(i)
            read(1,*) typat(i,j),Rx(i,j),Ry(i,j),Rz(i,j)
            Rx(i,j) = (Rx(i,j) / ax) + 0.5
            Ry(i,j) = (Ry(i,j) / ay) + 0.5
            Rz(i,j) = (Rz(i,j) / az) + 0.5            
          enddo
        enddo

 100    continue

        ! Find all dx,dy,dz values from the wrapped data
        do j=1,natoms(1)
          dx(1,j) = 0.0
          dy(1,j) = 0.0
          dz(1,j) = 0.0
        enddo
        do i=2,n
          do j=1,natoms(i)
            dRx = Rx(i,j) - Rx(i-1,j)
            dRy = Ry(i,j) - Ry(i-1,j)
            dRz = Rz(i,j) - Rz(i-1,j)
            if (dRx.ge.ax/2) then
                print *, "WARNING: dRx has exceeded a/2"
            endif
            if (dRy.ge.ay/2) then
                print *, "WARNING: dRy has exceeded a/2"
            endif
            if (dRz.ge.az/2) then
                print *, "WARNING: dRz has exceeded a/2"
            endif
            dx(i,j) = dx(i-1,j) + (dRx - anint(dRx))
            dy(i,j) = dy(i-1,j) + (dRy - anint(dRy))
            dz(i,j) = dz(i-1,j) + (dRz - anint(dRz))
          enddo
        enddo

!$omp  parallel do
!$omp& default(private)
!$omp& shared(natoms,timesteps)

        ! Unwrap the data
        do i=1,n
          if (natoms(i).eq.0) then
            goto 80
          endif
          do j=1,natoms(i)
            ! Calculate the unwrapped coordinate
            Rx(i,j) = Rx(1,j) + dx(i,j)
            Ry(i,j) = Ry(1,j) + dy(i,j)
            Rz(i,j) = Rz(1,j) + dz(i,j)
          enddo
        enddo

 80     continue

!$omp  end parallel do

        ! Write to file
        do i=1,n
          if (natoms(i).eq.0) then
            goto 70
          endif
          write(2,*) natoms(i)
          write(2,*) timesteps(i)
          do j=1,natoms(i)
            ! Convert back to Cartesian coordinates (re-shift origin)
            Rx(i,j) = (Rx(i,j) - 0.5) * ax
            Ry(i,j) = (Ry(i,j) - 0.5) * ay
            Rz(i,j) = (Rz(i,j) - 0.5) * az
            ! Write unwrapped coordinates to file
            write(2,*) typat(i,j), Rx(i,j), Ry(i,j), Rz(i,j)
          enddo
        enddo

 70     continue

        close(1)
        close(2)

 90     continue

        END
