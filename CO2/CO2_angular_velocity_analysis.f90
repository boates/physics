!       program CO2_anglular_velocity_analysis.f90
!***********************************************************
!       Calculates CO2 angular velocity, use this to then
!       determine equilibrium bond angles from peaked
!       distributions of angles occuring at velocity maxima
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, maxBins, i, j, k, f
        parameter (maxSteps=100000, maxAtoms=1024, maxBins=1000)
        integer nC, nSteps, bin, nbins, norm(maxAtoms), norm_all
        real*8 time(maxSteps), theta(maxAtoms,maxSteps), d_theta(maxAtoms,maxSteps)
        real*8 eq_angle, eq_angles(maxAtoms,maxBins), eq_angles_all(maxBins), tstep
        character*532 fmt

        ! Retrieve user input
        write(*,*) 'Number of carbon atoms:'
        read(*,*) nC

        ! Set nbins
        nbins = 200

        ! Open angle trajectories and output files
        open(1,file='CO2_tracked_angles.dat',status='old',ERR=90)
        open(2,file='angular_velocities_CO2.dat')
        open(3,file='equilibrium_angles_individual_CO2.hist')
        open(4,file='equilibrium_angles_all_CO2.hist')

        ! Read in angular trajectories
        do i=1,maxSteps
          read(1,*,END=110) time(i), (theta(j,i),j=1,nC)
        enddo

 110    continue

        ! Determine time step in ps
        tstep = time(2) - time(1)

        ! Record the number of timesteps
        nSteps = i - 1

        ! Calculate the angular velocities (differentiate, fwd-difference)
        do i=1,nC
          do j=1,nSteps-1
            d_theta(i,j) = (theta(i,j+1) - theta(i,j)) / (time(j+1) - time(j))
          enddo
        enddo

        ! Determine equilibrium angles from d_theta**2 peaks
        do i=1,nC
          do j=2,nSteps-1
            if (d_theta(i,j-1)**2.lt.d_theta(i,j)**2) then
              if (d_theta(i,j+1)**2.lt.d_theta(i,j)**2) then

                ! Because of forward difference the eq. angle is ~:
                eq_angle = ( theta(i,j) + theta(i,j+1) ) / 2.0

                ! Add to histograms
                bin = int(eq_angle*(nbins/180.0) + 0.5)
                eq_angles(i,bin) = eq_angles(i,bin) + 1
                eq_angles_all(bin) = eq_angles_all(bin) + 1
                norm(i) = norm(i) + 1
                norm_all = norm_all + 1

              endif
            endif
          enddo
        enddo

        ! Normalize the equilibrium distributions
        do i=1,nbins
          eq_angles_all(i) = eq_angles_all(i) / norm_all
          do j=1,nC
            eq_angles(j,i) = eq_angles(j,i) / norm(j)
          enddo
        enddo

        ! Create the format string
        fmt = '('
        f = 1
        do j=1,nC
          fmt = fmt(1:f) // 'f16.8,x,'
          f = f + 8
        enddo
        fmt = fmt(1:f) // 'f16.8)'

        ! Write angular trajectories to file
        do i=1,nSteps
          write(2,fmt) (i-1)*tstep, (d_theta(j,i),j=1,nC)
        enddo

        ! Write equilibrium angle distributions to file
        do i=1,nbins
          write(3,fmt) (i-1)*180.0/nbins, (eq_angles(j,i),j=1,nC)
          write(4,*) (i-1)*180.0/nbins, eq_angles_all(i)
        enddo

 100    continue

        close(1)
        close(2)
        close(3)
        close(4)

 90     continue

        END
