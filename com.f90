!       program com.f
!***********************************************************
!       Calculate the center of mass of a number of
!       identically massed particles, over a number
!       of timesteps (get the COM trajectory)
!***********************************************************
        implicit none

        integer n, i, j
        parameter (n=100000)
        integer natom, timestep
        real*8 Rx, Ry, Rz, sumRx, sumRy, sumRz
        character typat
        character*128 fin

        write(6,*) 'Name of unwrapped xyz file'
        read(5,*) fin

        open(1,file=fin,status='old',ERR=90)
        open(2,file='com.xyz')
        open(3,file='com.dat')

        do i=1,n

          read(1,*,END=100) natom
          read(1,*) timestep
          write(2,*) '1'
          write(2,*) timestep

          sumRx = 0
          sumRy = 0
          sumRz = 0

          do j=1,natom

            read(1,*) typat, Rx, Ry, Rz

            sumRx = sumRx + Rx
            sumRy = sumRy + Ry
            sumRz = sumRz + Rz

          enddo
          
          sumRx = sumRx / natom
          sumRy = sumRy / natom
          sumRz = sumRz / natom

          write(2,*) 'COM', sumRx, sumRy, sumRz

          write(3,*) timestep,(sumRx**2+sumRy**2+sumRz**2)**(0.5)

        enddo

 100    continue

        close(1)
        close(2)
        close(3)

 90     continue

        END
