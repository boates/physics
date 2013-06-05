!       program molecular_centers.f90
!***********************************************************
!       Create an xyz file of molecule-COMs
!       based on the information in a cbn file
!       ASSUME PURE MOLECULAR SYSTEM
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, tstep, tracker(maxSteps,maxAtoms), pbc_round(3)
        real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms), nn(maxAtoms)
        real*8 ax, ay, az, dx, dy, dz, distance, halfbox
        real*8 input_value(3)
        real Xcom, Ycom, Zcom
        character*16 typat, d1, d2, d3, d4, d5

        ! Open input and output files
        open(1,file='TRAJEC.cbn',status='old',ERR=90)
        open(2,file='CENTERS.xyz')

        ! Read cbn header info
        read(1,*,END=100) 
        read(1,*) d1, d2, d3, ax
        read(1,*) d1, d2, d3, ay
        read(1,*) d1, d2, d3, az
        read(1,*) d1, d2, d3, natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        halfbox = ( (ax/2)**2 + (ay/2)**2 + (az/2)**2 )**0.5

        do i=1,maxSteps

          ! Read data, and bump the nn index up by 1
          do j=1,natom
            read(1,*,END=100) typat, X(j), Y(j), Z(j), nn(j)
            nn(j) = nn(j) + 1
          enddo

          write(2,*) natom/2
          write(2,*) i

          k = 0
          do j=1,natom
            
            if (tracker(i,nn(j)).ne.-1) then

              k = k + 1

              ! Use the minimum image convention
              dx = X(j) - X(nn(j))
              dy = Y(j) - Y(nn(j))
              dz = Z(j) - Z(nn(j))
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
              Xcom = X(j) - dx/2.0
              Ycom = Y(j) - dy/2.0
              Zcom = Z(j) - dz/2.0

              ! Check that the atom to atom distance is reasonable
              distance = ( dx**2 + dy**2 + dz**2 )**0.5
              if (distance.gt.halfbox) then
                write(*,*) "WARNING --- PROBLEM WITH PBC, exiting..."
                write(*,*) X(j), X(nn(j))
                write(*,*) Y(j), Y(nn(j))
                write(*,*) Z(j), Z(nn(j))
                write(*,*) dx, dy, dz, distance
                write(*,*) ax, ay, az, halfbox, i, j
                goto 90
              endif

              write(2,*) 'C ', Xcom*0.529177, Ycom*0.529177, Zcom*0.529177

            endif

            tracker(i,j) = -1

          enddo

          ! Check that number of molecular centers = natom/2
          if (k.ne.natom/2) then
            write(*,*) "WARNING --- Number of molecular centers != natom/2, exiting..."
            write(*,*) i, k
            goto 100
          endif

        enddo

 100    continue

        close(1)

        close(2)

 90     continue

        END
