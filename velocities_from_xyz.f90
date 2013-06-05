!       program velocities_from_xyz.f90
!***********************************************************
!       Create an vxyz file of velocities from an xyz file
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, tstep
        real*8 Xp(maxAtoms), Yp(maxAtoms), Zp(maxAtoms)
        real*8 Xc(maxAtoms), Yc(maxAtoms), Zc(maxAtoms)
        real*8 dx, dy, dz, ax, ay, az
        real vx, vy, vz
        real*8 delta_t, velocity, temp(maxAtoms), temp_sum
        real*8 input_value(3), mass, kB
        integer pbc_round(3)
        character*8 typat
        character*128 fxyz

        ! Retrieve user input
        write(*,*) 'Name of xyz file'
        read(*,*) fxyz
        write(*,*) 'tstep in fs'
        read(*,*) delta_t
        write(*,*) 'ax, ay, az (angstrom)'
        read(*,*) ax,ay,az

        ! Open input and output files
        open(1,file=fxyz,status='old',ERR=90)
        open(2,file='VELOCITIES.vxyz')
        open(3,file='temp_calc.dat')
        open(4,file='track_atom.dat')

        ! Convert timestep size from AU to ps, set kB, and mass
        delta_t = delta_t / 1000.0
        kB = 0.83144727 ! in amu*A^2/ps^2/K
        mass = 14.007

        ! Read10 in first configuration
        read(1,*,END=100) natom
        read(1,*) tstep
        do j=1,natom
          read(1,*) typat, Xp(j), Yp(j), Zp(j)
        enddo

        ! Begin calculations
        do i=1,maxSteps

          temp_sum = 0.0

          read(1,*,END=100) natom
          read(1,*) tstep
          write(2,*) natom
          write(2,*) tstep

          do j=1,natom

            read(1,*) typat, Xc(j), Yc(j), Zc(j)

            ! Use the minimum image convention
            dx = Xc(j) - Xp(j)
            dy = Yc(j) - Yp(j)
            dz = Zc(j) - Zp(j)
            input_value(1) = dx/ax
            input_value(2) = dy/ay
            input_value(3) = dz/az
            pbc_round(1) = int(input_value(1))
            pbc_round(2) = int(input_value(2))
            pbc_round(3) = int(input_value(3))
            if (abs(input_value(1)-pbc_round(1)).ge.0.5) then
              if (input_value(1).gt.0) pbc_round(1) = pbc_round(1) + 1
              if (input_value(1).lt.0) pbc_round(1) = pbc_round(1) - 1
            endif
            if (abs(input_value(2)-pbc_round(2)).ge.0.5) then
              if (input_value(2).gt.0) pbc_round(2) = pbc_round(2) + 1
              if (input_value(2).lt.0) pbc_round(2) = pbc_round(2) - 1
            endif
            if (abs(input_value(3)-pbc_round(3)).ge.0.5) then
              if (input_value(3).gt.0) pbc_round(3) = pbc_round(3) + 1
              if (input_value(3).lt.0) pbc_round(3) = pbc_round(3) - 1
            endif
            dx = dx - ax*pbc_round(1)
            dy = dy - ay*pbc_round(2)
            dz = dz - az*pbc_round(3)

            ! Calculate velocities
            vx = dx / delta_t
            vy = dy / delta_t
            vz = dz / delta_t
            velocity = (dx**2 + dy**2 + dz**2)**0.5 / delta_t

            write(2,*) typat, vx, vy, vz

            if (j.eq.1) then
              write(4,*) i, vx, vy, vz, (vx**2+vy**2+vz**2)**0.5
            endif

            ! Calculate the temperature
            temp(j) = (mass * velocity**2) / (3.0 * kB)

            ! Rest previous coordinates as the current
            Xp(j) = Xc(j)
            Yp(j) = Yc(j)
            Zp(j) = Zc(j)

          enddo

          do j=1,natom
            temp_sum = temp_sum + temp(j)
          enddo

          write(3,*) temp_sum / natom
          
        enddo

 100    continue

        close(1)
        close(2)
        close(3)
        close(4)

 90     continue

        END

