!       program CO2_dist_angles.f90
!***********************************************************
!       Calculate a angular distribution function.
!       Output file is probability of finding other
!       molecules at a given angular "distance" away
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, maxBins, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024,maxBins=1000)
        integer nC, count, nBins, bin, norm
        real*8 time(maxSteps), angles(maxAtoms,maxSteps)
        real*8 da, angle_hist(maxBins)

        ! Retrieve user input
        write(*,*) 'Number of carbon atoms'
        read(*,*) nC

        ! Set the number of bins
        nBins = 360

        ! Open .cnn and output files
        open(1,file='CO2_tracked_angles.dat',status='old',ERR=90)
        open(2,file='ADF.dat')

        ! Perform analysis
        do i=1,maxSteps

          ! Read in angles for current timestep
          read(1,*,END=100) time(i), (angles(j,i),j=1,nC)
          count = 0

          ! Calculate angular difference for molecules
          do j=1,nC
            do k=j+1,nC
              if ((angles(j,i).ne.0.0).and.(angles(k,i).ne.0.0)) then
                da = angles(j,i) - angles(k,i)
                bin = int(da*nBins/180.0 + 0.5)
                angle_hist(bin) = angle_hist(bin) + 1
                norm = norm + 1
              endif
            enddo
          enddo

        enddo

 100    continue

        ! Write to output file
        do j=1,nBins
          angle_hist(j) = angle_hist(j) / norm
          write(2,'(f10.6,x,f10.6)') (j*180.0/nBins), angle_hist(j)
        enddo

        close(1)
        close(2)

        goto 80

 90     continue

        write(*,*) 'Make sure CO2_tracked_angles.dat file is present'

 80     continue

        END
