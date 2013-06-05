!       program CO2_molecule_detector.f90
!****************************************************************
!
!       Locates CO2 molecules from .cnn files and write to a
!       new file called a ".mol" file with cbn-ish form:
!
!       Output is in ANGSTROMS!
!
!       C x y z 1st_nearest_O 2nd_nearest_O
!       O x y z C_bonded_to other_O_bonded_to_same_C
!
!       Indicies of atoms start at ONE!
!
!       If atom not in defined molecule:
!       C x y z 0 0
!       O x y z 0 0
!
!****************************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=128)
        integer natom, nneighbors, nC, nO, mol(maxAtoms,maxAtoms)
        integer nn(maxAtoms,maxAtoms-1), nn1, nn2, c, o
        real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms), ax, ay, az
        character*8 typat(maxAtoms), d1, d2, d3
        character*128 fcnn

        ! Retrieve user input
        write(*,*) 'Name of TRAJEC.cnn file'
        read(*,*) fcnn
        write(*,*) 'Number of carbon atoms'
        read(*,*) nC
        write(*,*) 'Number of oxygen atoms'
        read(*,*) nO

        ! Open .cnn and output files
        open(1,file=fcnn,status='old',ERR=90)
        open(2,file='TRAJEC.mol')

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

        ! Make sure nC and nO consisent with natom
        if (natom.ne.nC+nO) then
          write(*,*) 'natom != nC + nO given --- exiting...'
          goto 95
        endif

        ! Write .mol header
        write(2,*) '# natom = ', natom
        write(2,*) '# nC =', nC
        write(2,*) '# nO =', nO
        write(2,*) '# units = angstroms'
        write(2,*) '# ax =', ax
        write(2,*) '# ay =', ay
        write(2,*) '# az =', az
        write(2,*) '# C x y z 1st_nearest_O 2nd_nearest_O'
        write(2,*) '# O x y z C_bonded_to other_O_bonded_to_same_C'
        write(2,*) '# If atom not a in molecule the last two columns are: 0 0'

        ! Perform analysis
        do i=1,maxSteps

          ! Read snapshot coordinates and neighbor lists
          do j=1,natom
            read(1,*,END=100) typat(j), X(j), Y(j), Z(j), (nn(j,k),k=1,nneighbors)
            X(j) = X(j)
            Y(j) = Y(j)
            Z(j) = Z(j)
          enddo

          do j=1,nC

            ! Make sure current atom is carbon
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

                  ! Make sure indices make sense for C and O's
                  if ((j.gt.nC).or.(nn1.le.nC).or.(nn2.le.nC)) then
                    write(*,*) 'C index, nC:', j, nC
                    write(*,*) 'O indicies, nO:', nn1, nn2, nO
                    write(*,*) 'At least one of these does not make sense --- exiting...'
                    goto 100
                  endif

                  mol(j,1) = j    ! Current carbon
                  mol(j,2) = nn1  ! Closest oxygen
                  mol(j,3) = nn2  ! 2nd closest oxygen

                endif
              endif
            endif
 210        continue
          enddo

          !========================================================!
          ! Entire current snapshot has been scanned for molecules !
          ! Molecule array (mol(i,j)) is ready to go at this point !
          ! mol(i,j) gives the indices of the carbon and 2 oxygens !
          ! in the i'th molecule                                   !
          !========================================================!

 5        format(A1,x,f12.8,x,f12.8,x,f12.8,x,x,i2,x,x,i2)

          ! Write carbon lines to .mol file
          do j=1,nC
            write(2,5) typat(j), X(mol(j,1)), Y(mol(j,1)), Z(mol(j,1)), mol(j,2), mol(j,3)
          enddo

          ! Write oxygen lines to .mol file
          do k=nC+1,natom
            do j=1,nC
              if (k.eq.mol(j,2)) then
                c = j
                o = mol(j,3)
              elseif (k.eq.mol(j,3)) then
                c = j
                o = mol(j,2)
              endif
            enddo
            write(2,5) typat(j), X(k), Y(k), Z(k), c, o
          enddo

        enddo

 100    continue

        close(2)

 95     continue

        close(1)

 90     continue

        END
