!       program vacf.f90
!***********************************************************
!       Calculate unaveraged VACF
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j
        parameter (maxSteps=10000,maxAtoms=1000)
        integer natom, tstep
        real*8 X(maxSteps,maxAtoms)
        real*8 Y(maxSteps,maxAtoms)
        real*8 Z(maxSteps,maxAtoms)
        real*8 dx, dy, dz, dt
        real*8 vx0(maxAtoms), vy0(maxAtoms), vz0(maxAtoms)
        real*8 vx, vy, vz, v0v0(maxAtoms), vacf(maxSteps)
        character*128 fin, typat

        ! Retrieve user input
        write(*,*) 'Name of unwrapped xyz file'
        read(*,*) fin
        write(*,*) 'timestep in a.u.'
        read(*,*) dt

        ! Convert a.u. to picoseconds
        dt = dt * 2.41880E-05 

        open(1,file=fin,status='old',ERR=90)
        open(2,file='vacf.dat')

        ! Read in configurations
        do i=1,maxSteps
          read(1,*,END=110) natom
          read(1,*) tstep
          do j=1,natom
            read(1,*) typat, X(i,j), Y(i,j), Z(i,j)
          enddo
          tstep = i
        enddo

 110    continue

        close(1)

        ! Calculate v(0) values
        do j=1,natom
          vx0(j)  = ( X(2,j) - X(1,j) ) / dt
          vy0(j)  = ( Y(2,j) - Y(1,j) ) / dt
          vz0(j)  = ( Z(2,j) - Z(1,j) ) / dt
          v0v0(j) = vx0(j)**2 + vy0(j)**2 + vz0(j)**2
        enddo

        ! Calculate the VACF for all configurations
        do i=2,tstep

          vacf(i-1) = 0.0

          do j=1,natom

            dx = X(i,j) - X(i-1,j)
            dy = Y(i,j) - Y(i-1,j)
            dz = Z(i,j) - Z(i-1,j)

            vx = dx / dt
            vy = dy / dt
            vz = dz / dt

            vacf(i-1) = vacf(i-1) + (vx*vx0(j)+vy*vy0(j)+vz*vz0(j)) / v0v0(j)

          enddo

          vacf(i-1) = vacf(i-1) / natom

        enddo

        ! Write VACF to file
        do i=1,tstep
          write(2,*) dt*(i-1), vacf(i)
        enddo

        close(2)

 90     continue

        END
