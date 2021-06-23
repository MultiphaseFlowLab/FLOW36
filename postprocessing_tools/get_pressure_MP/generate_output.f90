subroutine generate_output(nstep)

use commondata
use fields

integer :: nstep,nfields,numx,numy,numz
integer :: i,k,j

character(len=40) :: namefile
character(len=80) :: buffer
character(len=10) :: lf
character(len=8) :: str1,str2,str3
character(len=16) :: str4

! end of line character
lf=achar(10)

! fields included
nfields=2+phi_flag

numx=nx
numy=ny
numz=nz



write(namefile,'(a,i8.8,a)') './output/OUTPAR_',nstep,'.vtk'

 open(66,file=trim(namefile),status='new',form='unformatted',access='stream',convert='big_endian')

 ! start writing vtk file

 ! write header
 buffer='# vtk DataFile Version 3.0'//lf
 write(66) trim(buffer)
 buffer='Phase field'//lf
 write(66) trim(buffer)
 buffer='BINARY'//lf
 write(66) trim(buffer)
 buffer='DATASET RECTILINEAR_GRID'//lf
 write(66) trim(buffer)

 !write grid
 write(str1(1:8),'(i8)') numx
 write(str2(1:8),'(i8)') numy
 write(str3(1:8),'(i8)') numz
 buffer='DIMENSIONS '//str1//str2//str3//lf
 write(66) trim(buffer)
 buffer='X_COORDINATES '//str1//'  double'//lf ;
 write(66) trim(buffer)
 do i=1,nx
  write(66) x(i)
 enddo
 buffer='Y_COORDINATES '//str2//'  double'//lf ;
 write(66) trim(buffer)
 do j=1,ny
  write(66) y(j)
 enddo
 buffer='Z_COORDINATES '//str3//'  double'//lf ;
 write(66) trim(buffer)
 do k=1,nz
  write(66) z(k)
 enddo

 ! write content (data format)
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer='POINT_DATA '//str4//lf
 write(66) trim(buffer)
 write(str1(1:8),'(i8)') nfields
 buffer='FIELD FieldData '//str1//lf
 write(66) trim(buffer)

  ! write u field
  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'U 3 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=1,nz
   do j=1,ny
    do i=1,nx
     write(66) u(i,k,j),v(i,k,j),w(i,k,j)
    enddo
   enddo
  enddo

 ! write phi field
  if(phi_flag.eq.1)then
   write(str4(1:16),'(i16)') numx*numy*numz
   buffer = 'PHI 1 '//str4//' double'//lf
   write(66) trim(buffer)
   do k=1,nz
    do j=1,ny
     do i=1,nx
      write(66) phi(i,k,j)
     enddo
    enddo
   enddo
  endif

 ! write press field
  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'PRESS 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=1,nz
   do j=1,ny
    do i=1,nx
     write(66) press(i,k,j)
    enddo
   enddo
  enddo

 buffer=lf
 write(66) trim(buffer)

 close(66,status='keep')

return
end
