!       program CO2_orientation.f90
!****************************************************************
!       Calculates CO2 molecule releative orientations (part 1)
!****************************************************************
        implicit none
        
        integer maxSteps, maxAtoms, maxBins, i, j, k
        parameter (maxSteps=100000,maxAtoms=128, maxBins=1000)
        integer natom, nneighbors, nC, nCneighbors, Nbins, norm, count, bin
        integer Nnear, amount, Nbin, Mbin, OObin, MNnorm
        integer pbc_round(3), nn(maxAtoms,maxAtoms-1), nn1, nn2, nnMOL
        integer nnC(maxAtoms,maxAtoms), mol(maxAtoms,3)
        real*8 r, dx, dy, dz, ax, ay, az
        real*8 Amag, Bmag, Nmag, Mmag, kAmag, kBmag, kNmag, kMmag
        real*8 A(3), B(3), OO(3), M(3), N(3), kA(3), kB(3), kOO(3), kM(3), kN(3)
        real*8 Mangle, Nangle, Mangle_hist(maxBins), Nangle_hist(maxBins)
        real*8 OOangle, OOangle_hist(maxBins)
        real*8 angle, angle_hist(maxBins)
        real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
        real*8 com(maxAtoms,3), Cmass, Omass, xCOM, yCOM, zCOM
        real*8 input_value(3), pi, dC(maxAtoms,maxAtoms)
        character*8 typat(maxAtoms), d1, d2, d3, d4, d5
        character*128 fin, finC, label

        ! Retrieve user input
        write(*,*) 'Name of TRAJEC.cnn file'
        read(*,*) fin
        write(*,*) 'Name of TRAJEC_C.cnn file of at least same length'
        read(*,*) finC
        write(*,*) 'Calculate for how many nearest molecules? (1, 2, ..., 0=all)'
        read(*,*) Nnear
        write(*,*) 'Number of bins'
        read(*,*) Nbins

        write(label,'(I0)') Nnear
        ! Open .cnn and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file=finC,status='old',ERR=95)
        open(3,file='CO2_com_orientation_'//trim(label)//'.dat')
        open(4,file='CO2_plane_orientation_'//trim(label)//'.dat')
        open(5,file='CO2_molecule_fraction.dat')
!        open(6,file='CO2_molecule_angles.dat')
!        open(7,file='CO2_com.xyz')
!        open(8,file='oxygen_axes_angles.dat')

        ! Read in TRAJEC.cnn header
        read(1,*,END=100) 
        read(1,*) d1, d2, d3, ax
        read(1,*) d1, d2, d3, ay
        read(1,*) d1, d2, d3, az
        read(1,*) d1, d2, d3, natom
        read(1,*) d1, d2, d3, nneighbors
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        ! Read in TRAJEC_C.cnn header
        read(2,*,END=100) 
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) d1, d2, d3, nC
        read(2,*) d1, d2, d3, nCneighbors
        read(2,*) 
        read(2,*) 
        read(2,*) 
        read(2,*) 

        ! Define C and O atomic masses (in amu) (and pi)
        Cmass = 12.0107
        Omass = 15.9994
        pi = 3.141592653590

        ! Perform analysis
        do i=1,maxSteps

          ! Read snapshot coordinates and neighbor lists
          do j=1,natom
            read(1,*,END=110) typat(j), X(j), Y(j), Z(j), (nn(j,k),k=1,nneighbors)
          enddo

          ! Read carbon neighbor lists
          do j=1,nC
            read(2,*,END=110) d1, d2, d3, d4, (nnC(j,k),k=1,nCneighbors)
          enddo

          ! Restart molecule count for each timestep
          count = 0

          do j=1,natom

            ! If current atom is carbon
            if (typat(j).eq.'C') then
              nn1 = nn(j,1) + 1
              nn2 = nn(j,2) + 1
              ! If this carbon's two nearest neighbors are oxygen
              if (typat(nn1).eq.'O'.and.typat(nn2).eq.'O') then
                ! If this carbon is those two oxygen's 1st nearest neighbor
                if (nn(nn1,1).eq.j-1.and.nn(nn2,1).eq.j-1) then
                  ! Check to make sure these oxygens aren't one of another carbon's
                  ! two nearest neighbors (avoid extended structures, O-C-O-C-O... etc.)
                  do k=1,natom
                    if (typat(k).eq.'C') then
                      if (j.eq.k) then
                        goto 200
                      endif
                      if (nn1.eq.nn(k,1)+1.or.nn1.eq.nn(k,2)+1) then
                        goto 210
                      endif
                      if (nn2.eq.nn(k,1)+1.or.nn2.eq.nn(k,2)+1) then
                        goto 210
                      endif
                    endif
                    200 continue
                  enddo

                  ! count the number of valid molecules found for current timestep
                  count = count + 1

                  !===============================================================!
                  ! At this point, we should have found a carbon with first two   !
                  ! neighbors being oxygen and which the oxygens are no other     !
                  ! carbon's two nearest neighbors, but are both have the current !
                  ! carbon as their first nearest neighbor. Let's continue.       !
                  !===============================================================!

                  ! Append the molecule to the molecule array
                  mol(count,1) = j    ! Current carbon
                  mol(count,2) = nn1  ! Closest oxygen
                  mol(count,3) = nn2  ! 2nd closest oxygen

                  ! Make sure indices make sense for C and O's
                  if (j.gt.nC.or.nn1.le.nC.or.nn2.le.nC) then
                    write(*,*) 'C index:', j,nC
                    write(*,*) 'O indicies:', nn1, nn2
!                    write(*,*) 'At least one of these does not make sense, EXITING...'
!                    goto 100
                  endif

                endif
              endif
            endif

          210 continue
          enddo

          !========================================================!
          ! Entire current snapshot has been scanned for molecules !
          ! Molecule array (mol(i,j)) is ready to go at this point !
          ! mol(i,j) gives the indices of the carbon and 2 oxygens !
          ! in the i'th molecule                                   !
          !========================================================!

          ! Check for too many CO2 molecules
          if (count.gt.nC) then
            write(*,*) 'WARNING: Number of CO2 molecules > detected nC --- EXITING...'
            write(*,*) 'i.e.',count,'>',nC
            goto 100
          endif

          ! Calculate center of masses of the molecules                                                                           
          ! com(j,i) gives the i'th component of the j'th molecule's position (mol(j,k))                                          
!          write(7,*) count
!          write(7,*) i
!          do j=1,count
!            xCOM = Cmass * X(mol(j,1)) + Omass * ( X(mol(j,2)) + X(mol(j,3)) )
!            yCOM = Cmass * Y(mol(j,1)) + Omass * ( Y(mol(j,2)) + Y(mol(j,3)) )
!            zCOM = Cmass * Z(mol(j,1)) + Omass * ( Z(mol(j,2)) + Z(mol(j,3)) )
!            com(j,1) = xCOM / (Cmass * 2*Omass)
!            com(j,2) = yCOM / (Cmass * 2*Omass)
!            com(j,3) = zCOM / (Cmass * 2*Omass)
!            write(7,*) 'C', com(j,1), com(j,2), com(j,3)
!          enddo

          ! Output fraction of CO2 molecules to carbon atoms
          write(5,*) i, real(count) / real(nC)

          ! Check for amount of molecules to calculate out to
          if (Nnear.gt.count.or.Nnear.eq.0) then
            amount = count - 1
          else
            amount = Nnear
          endif

          ! Loop over all molecules
          do j=1,count
            ! Calculate self vectors for j'th molecule
            !=========== Use the minimum image convention ===================!
            dx = X(mol(j,2)) - X(mol(j,1))
            dy = Y(mol(j,2)) - Y(mol(j,1))
            dz = Z(mol(j,2)) - Z(mol(j,1))
            input_value(1) = dx / ax
            input_value(2) = dy / ay
            input_value(3) = dz / az
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
            A(1) = dx - ax*pbc_round(1)
            A(2) = dy - ay*pbc_round(2)
            A(3) = dz - az*pbc_round(3)
       
            ! Final result for distance between CO's        
            r = ( A(1)**2 + A(2)**2 + A(3)**2 )**0.5

            ! Check that the distance makes sense
            if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
              write(*,*) "WARNING: Problem with PBC --- EXITING..."
              goto 100
            endif
            !================================================================!
            dx = X(mol(j,3)) - X(mol(j,1))
            dy = Y(mol(j,3)) - Y(mol(j,1))
            dz = Z(mol(j,3)) - Z(mol(j,1))
            input_value(1) = dx / ax
            input_value(2) = dy / ay
            input_value(3) = dz / az
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
            B(1) = dx - ax*pbc_round(1)
            B(2) = dy - ay*pbc_round(2)
            B(3) = dz - az*pbc_round(3)
       
            ! Final result for distance between CO's        
            r = ( B(1)**2 + B(2)**2 + B(3)**2 )**0.5

            ! Check that the distance makes sense
            if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
              write(*,*) "WARNING: Problem with PBC --- EXITING..."
              goto 100
            endif
            !================================================================!
            dx = X(mol(j,3)) - X(mol(j,2))
            dy = Y(mol(j,3)) - Y(mol(j,2))
            dz = Z(mol(j,3)) - Z(mol(j,2))
            input_value(1) = dx / ax
            input_value(2) = dy / ay
            input_value(3) = dz / az
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
            OO(1) = dx - ax*pbc_round(1)
            OO(2) = dy - ay*pbc_round(2)
            OO(3) = dz - az*pbc_round(3)
       
            ! Final result for distance between O's        
            r = ( OO(1)**2 + OO(2)**2 + OO(3)**2 )**0.5

            ! Check that the distance makes sense
            if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
              write(*,*) "WARNING: Problem with PBC --- EXITING..."
              goto 100
            endif
            !================================================================!

            ! Calculate central axis M(i) and plane normal N(i) vectors for j'th molecule
            Mmag = ((A(1)+B(1))**2+(A(2)+B(2))**2+(A(3)+B(3))**2)**0.5
            M(1) = (A(1)+B(1)) / Mmag
            M(2) = (A(2)+B(2)) / Mmag
            M(3) = (A(3)+B(3)) / Mmag

            N(1) = A(2)*B(3) - A(3)*B(2)
            N(2) = A(3)*B(1) - A(1)*B(3)
            N(3) = A(1)*B(2) - A(2)*B(1)
            Nmag = ( N(1)**2 + N(2)**2 + N(3)**2 )**0.5
            N(1) = N(1) / Nmag
            N(2) = N(2) / Nmag
            N(3) = N(3) / Nmag

            ! Calculate CO2 internal angles
            Amag = ( A(1)**2 + A(2)**2 + A(3)**2 )**0.5
            Bmag = ( B(1)**2 + B(2)**2 + B(3)**2 )**0.5
            angle = acos( ( A(1)*B(1) + A(2)*B(2) + A(3)*B(3) ) / ( Amag*Bmag ) )
            angle = angle * (180.0/3.141592653590)

            bin = int(angle*Nbins/180.0 + 0.5)
            angle_hist(bin) = angle_hist(bin) + 1
            norm = norm + 1

            ! Calculate out to requested amount of molecules
            do k=1,amount
              ! k'th nearest carbon to j'th molecule's carbon
              nnMOL = nnC(mol(j,1),k) + 1
              ! Calculate self vectors for nnMOL
              !=========== Use the minimum image convention ===================!
              dx = X(mol(nnMOL,2)) - X(nnMOL)
              dy = Y(mol(nnMOL,2)) - Y(nnMOL)
              dz = Z(mol(nnMOL,2)) - Z(nnMOL)
              input_value(1) = dx / ax
              input_value(2) = dy / ay
              input_value(3) = dz / az
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
              kA(1) = dx - ax*pbc_round(1)
              kA(2) = dy - ay*pbc_round(2)
              kA(3) = dz - az*pbc_round(3)
       
              ! Final result for distance between CO's        
              r = ( kA(1)**2 + kA(2)**2 + kA(3)**2 )**0.5

              ! Check that the distance makes sense
              if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
                write(*,*) "WARNING: Problem with PBC --- EXITING..."
                goto 100
              endif
              !================================================================!
              dx = X(mol(nnMOL,3)) - X(nnMOL)
              dy = Y(mol(nnMOL,3)) - Y(nnMOL)
              dz = Z(mol(nnMOL,3)) - Z(nnMOL)
              input_value(1) = dx / ax
              input_value(2) = dy / ay
              input_value(3) = dz / az
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
              kB(1) = dx - ax*pbc_round(1)
              kB(2) = dy - ay*pbc_round(2)
              kB(3) = dz - az*pbc_round(3)
       
              ! Final result for distance between CO's        
              r = ( kB(1)**2 + kB(2)**2 + kB(3)**2 )**0.5

              ! Check that the distance makes sense
              if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
                write(*,*) "WARNING: Problem with PBC --- EXITING..."
                goto 100
              endif
              !================================================================!
              dx = X(mol(nnMOL,3)) - X(mol(nnMOL,2))
              dy = Y(mol(nnMOL,3)) - Y(mol(nnMOL,2))
              dz = Z(mol(nnMOL,3)) - Z(mol(nnMOL,2))
              input_value(1) = dx / ax
              input_value(2) = dy / ay
              input_value(3) = dz / az
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
              kOO(1) = dx - ax*pbc_round(1)
              kOO(2) = dy - ay*pbc_round(2)
              kOO(3) = dz - az*pbc_round(3)
       
              ! Final result for distance between O's        
              r = ( kOO(1)**2 + kOO(2)**2 + kOO(3)**2 )**0.5

              ! Check that the distance makes sense
              if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
                write(*,*) "WARNING: Problem with PBC --- EXITING..."
                goto 100
              endif
              !================================================================!

              ! Now have both C-O vectors for j'th molecule and its k'th nearest molecule
              
              ! Calculate central axis M(i) and plane normal N(i) vectors for j'th molecule's
              ! k'th nearest molecule
              kMmag = ((kA(1)+kB(1))**2+(kA(2)+kB(2))**2+(kA(3)+kB(3))**2)**0.5
              kM(1) = (kA(1)+kB(1)) / kMmag
              kM(2) = (kA(2)+kB(2)) / kMmag
              kM(3) = (kA(3)+kB(3)) / kMmag

              kN(1) = kA(2)*kB(3) - kA(3)*kB(2)
              kN(2) = kA(3)*kB(1) - kA(1)*kB(3)
              kN(3) = kA(1)*kB(2) - kA(2)*kB(1)
              kNmag = ( kN(1)**2 + kN(2)**2 + kN(3)**2 )**0.5
              kN(1) = kN(1) / kNmag
              kN(2) = kN(2) / kNmag
              kN(3) = kN(3) / kNmag

              ! Caclulate angles between M and N vectors
              ! Note: M and N are already both normalized to unity magnitude
              Mangle = acos( M(1)*kM(1) + M(2)*kM(2) + M(3)*kM(3) )
              Nangle = acos( N(1)*kN(1) + N(2)*kN(2) + N(3)*kN(3) )
!              OOangle = acos( OO(1)*kOO(1) + OO(2)*kOO(2) + OO(3)*kOO(3) )

              Mangle = Mangle * (180.0 / pi)
              Nangle = Nangle * (180.0 / pi)
!              OOangle = OOangle * (180.0 / pi)

              ! Store angle data in histograms
              Mbin = int(Mangle*Nbins/180.0 + 0.5)
              Mangle_hist(Mbin+1) = Mangle_hist(Mbin+1) + 1
              Nbin = int(Nangle*Nbins/180.0 + 0.5)
              Nangle_hist(Nbin+1) = Nangle_hist(Nbin+1) + 1
!              OObin = int(OOangle*Nbins/180.0 + 0.5)
!              OOangle_hist(OObin+1) = OOangle_hist(OObin+1) + 1

              MNnorm = MNnorm + 1

            enddo

          enddo
        enddo

 110    continue

        do i=1,Nbins+1
          Mangle_hist(i) = Mangle_hist(i) / MNnorm
          Nangle_hist(i) = Nangle_hist(i) / MNnorm
!          OOangle_hist(i) = OOangle_hist(i) / MNnorm
!          angle_hist(i)  = angle_hist(i)  / norm
        enddo

        do i=1,Nbins+1          
          write(3,*) ((i-1)*180.0/Nbins), Mangle_hist(i)
          write(4,*) ((i-1)*180.0/Nbins), Nangle_hist(i)
!          write(6,*) ((i-1)*180.0/Nbins), angle_hist(i)
!          write(8,*) ((i-1)*180.0/Nbins), OOangle_hist(i)
        enddo

        close(3)
        close(4)
        close(5)
!        close(6)
!        close(7)
!        close(8)

 100    continue

        close(2)

 95     continue

        close(1)

 90     continue

        END
