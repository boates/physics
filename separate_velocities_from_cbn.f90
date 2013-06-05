!       program separate_velocities_from_cbn.f90
!***********************************************************
!       Create an vxyz files of velocities from an cbn file
!       One for bonded, one for non-bonded
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, tstep
        real*8 Xp(maxAtoms), Yp(maxAtoms), Zp(maxAtoms), NNp(maxAtoms)
        real*8 Xc(maxAtoms), Yc(maxAtoms), Zc(maxAtoms), NNc(maxAtoms)
        real*8 dx, dy, dz, ax, ay, az
        real vx, vy, vz
        real*8 delta_t, velocity
        real*8 temp(maxAtoms), temp_B(maxAtoms), temp_NB(maxAtoms)
        real*8 temp_sum, temp_sum_B, temp_sum_NB
        real*8 input_value(3), mass, kB
        integer pbc_round(3), Nchain, Nbonded
        character*8 typat
        character*128 fcbn, d1, d2, d3

        ! Retrieve user input
        write(*,*) 'Name of cbn file'
        read(*,*) fcbn
        write(*,*) 'tstep in AU'
        read(*,*) delta_t

        ! Open input and output files
        open(1,file=fcbn,status='old',ERR=90)
        open(2,file='B_temp_calc.dat')
        open(3,file='NB_temp_calc.dat')
        open(4,file='diff_temp_calc.dat')
        open(5,file='temp_calc.dat')

        ! Read cbn header
        read(1,*) 
        read(1,*) d1, d2, d3, ax
        read(1,*) d1, d2, d3, ay
        read(1,*) d1, d2, d3, az
        read(1,*) d1, d2, d3, natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        ! Convert lattice parameters to angtroms
        ax = ax*0.529177
        ay = ay*0.529177
        az = az*0.529177

        ! Convert timestep size from AU to ps, set kB, and mass
        delta_t = delta_t * 2.4188E-17
        delta_t = delta_t / 1.0E-12
        kB = 0.83144727 ! in amu*A^2/ps^2/K
        mass = 14.007

        ! Read in first configuration
        do j=1,natom
          read(1,*,END=100) typat, Xp(j), Yp(j), Zp(j), NNp(j)
          Xp(j) = Xp(j)*0.529177
          Yp(j) = Yp(j)*0.529177
          Zp(j) = Zp(j)*0.529177
        enddo

        ! Begin calculations
        do i=1,maxSteps

          Nbonded = 0
          Nchain = 0
          temp_sum = 0.0
          temp_sum_B = 0.0
          temp_sum_NB = 0.0

          do j=1,natom

            read(1,*,END=100) typat, Xc(j), Yc(j), Zc(j), NNc(j)

            Xc(j) = Xc(j)*0.529177
            Yc(j) = Yc(j)*0.529177
            Zc(j) = Zc(j)*0.529177

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

            ! Calculate the temperature
            temp(j) = (mass * velocity**2) / (3.0 * kB)
            if (NNc(j).ne.-1) then
              temp_B(j) = (mass * velocity**2) / (3.0 * kB)
              Nbonded = Nbonded + 1
            endif
            if (NNc(j).eq.-1) then
              temp_NB(j) = (mass * velocity**2) / (3.0 * kB)
              Nchain = Nchain + 1
            endif

            ! Rest previous coordinates as the current
            Xp(j) = Xc(j)
            Yp(j) = Yc(j)
            Zp(j) = Zc(j)
            NNp(j) = NNc(j)

          enddo

          if (Nbonded+Nchain.ne.natom) then
            write(*,*) 'WARNING: Nbonded + Nchain != natom at step', i+1
            write(*,*) 'Nbonded = ', Nbonded
            write(*,*) 'Nchain = ', Nchain
            write(*,*) 'natom = ', natom
          endif

          do j=1,natom
            temp_sum = temp_sum + temp(j)
            temp_sum_B = temp_sum_B + temp_B(j)
            temp_sum_NB = temp_sum_NB + temp_NB(j)
          enddo

          write(2,*) temp_sum_B / Nbonded
          write(3,*) temp_sum_NB / Nchain
          write(4,*) (temp_sum_NB / Nchain) - (temp_sum_B / Nbonded)
          write(5,*) temp_sum / natom
          
        enddo

 100    continue

        close(1)
        close(2)
        close(3)
        close(4)
        close(5)

 90     continue

        END

