!       program enthalpy.f90
!***********************************************************
!       Calculate the enthalpy from a pressure.dat,
!       energy.dat, and volume.dat
!***********************************************************
        implicit none
        
        integer maxSteps, i, j, k
        parameter (maxSteps=1000000)
        real*8 E(maxSteps), P(maxSteps), V, H(maxSteps)

        ! Open Required input and output files
        open(1,file='energy.dat',status='old',ERR=90)
        open(2,file='pressure.dat',status='old',ERR=90)
        open(3,file='volume.dat',status='old',ERR=90)
        read(2,*)
        read(3,*) V
        open(9,file='enthalpy.dat')

        ! Read in data an calculate the enthalpy in eV
        do i=1,maxSteps

          read(1,*,END=100) E(i)
          read(2,*,END=100) k, P(i)

          ! Convert from GPa to eV/angstrom^3
          P(i) = P(i) * 0.0062415097

          ! Calculate the enthalpy (eV)
          H(i) = E(i) + P(i)*V

          ! Write the enthalpy to file
          write(9,*) H(i)

        enddo

 100    continue

        close(1)
        close(2)
        close(3)
        close(9)

 90     continue

        END
