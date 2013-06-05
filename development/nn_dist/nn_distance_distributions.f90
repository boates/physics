!       program nn_distance_distributions.f
!***********************************************************
!       Translation of Isaac's nn_distance_distributions.py
!       into fortran for speed
!       Still pretty slow though...
!***********************************************************
        implicit none
        
        integer maxSteps, maxBins, maxAtoms
        parameter (maxSteps=1000000, maxBins=10000, maxAtoms=10000)
        integer i, j, s, p, o
        integer natom, Nneighbors, tstep, step_initial, Nsteps
        integer Nbins, stride, bin, count
        integer pbc_round(3), histogram(maxAtoms,maxBins)
        real*8 alat, ax, ay, az, bin_size
        real*8 r(3), dx, dy, dz, dr
        real*4 POS_array(maxAtoms,3), distance_array(maxAtoms,maxAtoms)
        real*8 input_value(3), distance, min_distance, max_distance
        real*8 min_distances(maxAtoms), min_distances_sort(maxAtoms)
        character*2 typat
        character*128 fxyz

        ! Get the xyz filename and alat from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fxyz
        write(*,*) 'Lattice constant (bohr)'
        read(*,*) alat
        write(*,*) 'Number of bins'
        read(*,*) Nbins
!        write(*,*) 'Snapshot stride'
!        read(*,*) stride

        ! Set and check the stride variable
        stride = 1
        if (stride.lt.1) then
          stride = 1
        endif

        ! Set the individual lattice constants
        ax = alat
        ay = alat
        az = alat

        ! Open the input and output files
        open(1,file=fxyz,status='old',ERR=90)
        open(2,file='nn_average.hist')
        write(2,*) '# bin (Angst), nn1, nn2, ...'

        ! Read in natom and tsteps from the xyz file
        do i=1,maxSteps
          read(1,*,END=110) natom
          read(1,*) tstep
          if (i.eq.1) then
             step_initial = tstep
          endif
          do j=1,natom
            read(1,*,END=110)
          enddo
        enddo

 110    continue

        close(1)

        ! Check the number of atoms
        if (natom.gt.maxAtoms) then
          write(*,*) 'natom =',natom,' exceeds maxAtoms =',maxAtoms
          write(*,*) 'Exiting...'
          goto 90
        endif

        ! Caclulate number of steps and neighbors
        Nsteps = tstep - step_initial + 1
        Nneighbors = natom - 1

        ! Histogram stuff
        min_distance = 0
        max_distance = alat / 2.0
        bin_size = (max_distance - min_distance) / Nbins

        ! Reopen the file
        open(1,file=fxyz,status='old',ERR=90)
        
        s = 0
        do while (s.lt.maxSteps)

          if (mod(s,100).eq.0.and.s.ne.0) then
            write(*,*) s*100.0/Nsteps,'%'
          endif

          read(1,*,END=100) natom
          read(1,*) tstep
          write(*,*) 'tstep=',tstep

          do p=1,natom
            read(1,*) typat, r(1), r(2), r(3)
            POS_array(p,1) = r(1)/0.5291772
            POS_array(p,2) = r(2)/0.5291772
            POS_array(p,3) = r(3)/0.5291772         
          enddo

          count = 1

          do p=1,natom

            do o=1,natom
              ! Calculate differences
              dx = POS_array(p,1) - POS_array(o,1)
              dy = POS_array(p,2) - POS_array(o,2)
              dz = POS_array(p,3) - POS_array(o,3)

              ! Use the minimum image convention
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

              distance = ( dx**2 + dy**2 + dz**2 )**0.5

              min_distances(o) = distance
            enddo

            ! Sort the minimum distances array to a new array
            do i=1,natom
              min_distances_sort(i) = minval(min_distances)
              min_distances(minloc(min_distances)) = maxval(min_distances)+1
            enddo

            do i=1,Nneighbors
              distance_array(count,i) = min_distances_sort(i+1)
            enddo
            count = count + 1
          enddo

          do i=1,count
            do j=1,Nneighbors
              if (distance_array(i,j).lt.max_distance) then
                bin = int( (distance_array(i,j) - min_distance)/bin_size )
                histogram(j,bin) = histogram(j,bin) + 1
              endif
            enddo
          enddo

          s = s + stride

        enddo

 100    continue

        close(1)

        write(*,*) 'natom, Nsteps, natom*Nsteps',natom,Nsteps,natom*Nsteps
!        namelist/output/tstep,time,P,V,T,E,H,msd

        ! Write the output file
        do i=1,Nbins
          write(2,*) (i*(bin_size/2.0)+min_distance)*0.529177,' ',histogram(j,:Nneighbors)/(natom*Nsteps)
!          do j=1,Nneighbors
!            write(2,*) histogram(j,i)/(natom*tstep),' ',
!          enddo
 !         write(2,*) ''
        enddo

        close(2)

 90     continue

        END
