! program nn_angles.f90
!****************************************************
! Calculate the angle between atom i and its j'th
! and k'th nearest neighbours (centered on i)
!
! Specify which type of atoms to center the
! calculation on (i.e. only on carbons a.k.a.
! O-C-O angles without C-O-O, C-C-O etc.)
!****************************************************

implicit none

integer maxSteps, maxAtoms, maxBins, i, j, k, count
parameter (maxSteps=100000,maxAtoms=1024,maxBins=1000)
integer natom, nnear, inner, outer
integer nn(maxAtoms,maxAtoms-1), nn2(maxAtoms,maxAtoms-1)
integer Nbins, bin, norm, pbc_round(3)
real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
real*8 ax, ay, az, dx, dy, dz, r
real*8 A(3), B(3), Amag, Bmag, dot, angle
real*8 input_value(3), angle_hist(maxBins)
character*2 typat(maxAtoms), all, s1, s2
character*128 fcnn, label1, label2, d1, d2, d3

! Retrieve user input
write(*,*) "TRAJEC.cnn file"
read(*,*) fcnn
write(*,*) "inner neighbour, outer neighbour"
read(*,*) inner, outer
write(*,*) "Number of bins"
read(*,*) Nbins
write(*,*) "All species? (y/n)"
read(*,*) all
if (all.eq.'n') then
  write(*,*) "Species 1, Species 2 (i.e. C, O for only O-C-O angles)"
  read(*,*) s1, s2
else
  s1 = 'A'
  s2 = 'A'
endif

write(label1,'(I0)') inner
write(label2,'(I0)') outer

! Open necessary files
open(1,file=fcnn,status='old',ERR=90)
open(2,file='angles_'//trim(s1)//trim(s2)//'_'//trim(label1)//'-'//trim(label2)//'.dat',status='unknown')

! Read the header
read(1,*)
read(1,*) d1, d2, d3, ax
read(1,*) d1, d2, d3, ay
read(1,*) d1, d2, d3, az
read(1,*) d1, d2, d3, natom
read(1,*) d1, d2, d3, nnear
read(1,*)
read(1,*)
read(1,*)
read(1,*)

! Loop over all timesteps in file
do i = 1, maxSteps

  ! Read in trajectories and neighbour list
  do j = 1, natom
    read(1,*,END=110) typat(j), X(j), Y(j), Z(j), (nn(j,k),k=1,nnear)
    if (all.ne.'n') then
      typat(j) = 'A'
    endif
    ! Shift the indices of the neighbour list to be fortran friendly
    do k = 1, nnear
      nn(j,k) = nn(j,k) + 1
    enddo
  enddo

  ! Calculate the angles for the desired species
  do j = 1, natom

    ! If the current atom is of species 1
    if (typat(j).eq.s1) then

      ! Create a "species 2 only" neighbour list
      count = 0
      do k = 1, nnear
        ! If a neighbour is species 2 add it to the new list
        if (typat( nn(j,k) ).eq.s2) then
          count = count + 1
          nn2(j,count) = nn(j,k)
        endif
      enddo

      ! Calculate the vector to inner
      dx = X(nn2(j,inner)) - X(j)
      dy = Y(nn2(j,inner)) - Y(j)
      dz = Z(nn2(j,inner)) - Z(j)
      input_value(1) = dx/ax
      input_value(2) = dy/ay
      input_value(3) = dz/az
      do k = 1, 3
        pbc_round(k) = int(input_value(k))
        if (abs(input_value(k)-pbc_round(k)).ge.0.5) then
          if (input_value(k).gt.0) pbc_round(k) = pbc_round(k) + 1
          if (input_value(k).lt.0) pbc_round(k) = pbc_round(k) - 1
        endif
      enddo
      A(1) = dx - ax*pbc_round(1)
      A(2) = dy - ay*pbc_round(2)
      A(3) = dz - az*pbc_round(3)
      r = ( A(1)**2 + A(2)**2 + A(3)**2 )**0.5  
      if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
        write(*,*) "WARNING --- PROBLEM WITH PBC"
        goto 100
      endif

      ! Calculate the vector to outer
      dx = X(nn2(j,outer)) - X(j)
      dy = Y(nn2(j,outer)) - Y(j)
      dz = Z(nn2(j,outer)) - Z(j)
      input_value(1) = dx/ax
      input_value(2) = dy/ay
      input_value(3) = dz/az
      do k = 1, 3
        pbc_round(k) = int(input_value(k))
        if (abs(input_value(k)-pbc_round(k)).ge.0.5) then
          if (input_value(k).gt.0) pbc_round(k) = pbc_round(k) + 1
          if (input_value(k).lt.0) pbc_round(k) = pbc_round(k) - 1
        endif
      enddo
      B(1) = dx - ax*pbc_round(1)
      B(2) = dy - ay*pbc_round(2)
      B(3) = dz - az*pbc_round(3)
      r = ( B(1)**2 + B(2)**2 + B(3)**2 )**0.5  
      if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
        write(*,*) "WARNING --- PROBLEM WITH PBC"
        goto 100
      endif

      ! Calculate the angle
      dot   = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      Amag  = ( A(1)**2 + A(2)**2 + A(3)**2 )**0.5
      Bmag  = ( B(1)**2 + B(2)**2 + B(3)**2 )**0.5
      angle = acos( dot / (Amag*Bmag) ) * (180.0/3.141592653590)
      bin   = int(angle*Nbins/180.0 + 0.5)
      angle_hist(bin) = angle_hist(bin) + 1
      norm  = norm + 1

    endif
      
  enddo
enddo

110 continue

! Write angles histogram to file
do i = 1, Nbins
  write(2,*) (i*180.0/Nbins), angle_hist(i) !/ norm
enddo

100 continue

close(1)
close(2)

90 continue

END
