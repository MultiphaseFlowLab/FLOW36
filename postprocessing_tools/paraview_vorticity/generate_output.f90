subroutine generate_output(nstep)

use commondata

implicit none

integer :: nstep,nfields,numx,numy,numz
integer :: i,j,k

character(len=40) :: namefile
character(len=80) :: buffer
character(len=10) :: lf
character(len=8) :: str1,str2,str3
character(len=16) :: str4

double precision, dimension(nx,nz,ny) :: up,vp,wp,omx,omy,omz,strx,stry,strz
double precision, dimension(nx/2+1,nz,ny,2) :: tmpc
double precision :: um,vm,wm


! end of line character
lf=achar(10)

! u,v,w fields + u',v',w' + omega x,y,z  + phi (if included)
nfields=3+3*flucflag+3*vorflag+3*strflag+phiflag+psiflag+tempflag


if(flucflag.eq.1)then
  ! calculate velocity fluctuations
  do k=1,nz
    um=0.0d0
    vm=0.0d0
    wm=0.0d0
    do j=1,ny
      do i=1,nx
        um=um+u(i,k,j)
        vm=vm+v(i,k,j)
        wm=wm+w(i,k,j)
      enddo
    enddo
    um=um/dble(nx*ny)
    vm=vm/dble(nx*ny)
    wm=wm/dble(nx*ny)
    up(:,k,:)=u(:,k,:)-um
    vp(:,k,:)=v(:,k,:)-vm
    wp(:,k,:)=w(:,k,:)-wm
  enddo
endif


if(vorflag.eq.1)then
  ! calculate vorticity
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)

  ! om x
  call dz(-vc,tmpc)
  do j=1,ny
    do i=1,nx/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)-ky(j)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)+ky(j)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,omx,0)

  ! om y
  call dz(uc,tmpc)
  do j=1,ny
    do i=1,nx/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)+kx(i)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)-kx(i)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,omy,0)

  ! om z
  do j=1,ny
    do i=1,nx/2+1
      tmpc(i,:,j,1)=-kx(i)*vc(i,:,j,2)+ky(j)*uc(i,:,j,2)
      tmpc(i,:,j,2)=+kx(i)*vc(i,:,j,1)-ky(j)*uc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,omz,0)

  ! normalize to wall units
  omx=omx/re
  omy=omy/re
  omz=omz/re
endif


if(strflag.eq.1) then
  ! calculate strain rate
  call phys_to_spectral(u,uc,0)
  call phys_to_spectral(v,vc,0)
  call phys_to_spectral(w,wc,0)

  ! str x
  call dz(vc,tmpc)
  do j=1,ny
    do i=1,nx/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)-ky(j)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)+ky(j)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,strx,0)
  strx=strx*0.5d0

  ! str y
  call dz(uc,tmpc)
  do j=1,ny
    do i=1,nx/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)-kx(i)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)+kx(i)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,stry,0)
  stry=stry*0.5d0

  ! str z
  do j=1,ny
    do i=1,nx/2+1
      tmpc(i,:,j,1)=-ky(j)*uc(i,:,j,2)-kx(i)*vc(i,:,j,2)
      tmpc(i,:,j,2)=+ky(j)*uc(i,:,j,1)+kx(i)*vc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,strz,0)
  strz=strz*0.5d0

  ! normalize to wu
  strx=strx/re
  stry=stry/re
  strz=strz/re
endif


numx=0
numy=0
numz=0
do i=x_start,x_end,dnx
 numx=numx+1
enddo
do j=y_start,y_end,dny
 numy=numy+1
enddo
do k=z_start,z_end,dnz
 numz=numz+1
enddo



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
 do i=x_start,x_end,dnx
  write(66) x(i)
 enddo
 buffer='Y_COORDINATES '//str2//'  double'//lf ;
 write(66) trim(buffer)
 do j=y_start,y_end,dny
  write(66) y(j)
 enddo
 buffer='Z_COORDINATES '//str3//'  double'//lf ;
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
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
 buffer = 'U 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) u(i,k,j)
   enddo
  enddo
 enddo

 ! write v field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'V 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) v(i,k,j)
   enddo
  enddo
 enddo

 ! write w field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'W 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) w(i,k,j)
   enddo
  enddo
 enddo


if(flucflag.eq.1)then
 ! write u' field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'UP 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) up(i,k,j)
   enddo
  enddo
 enddo

 ! write v' field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'VP 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) vp(i,k,j)
   enddo
  enddo
 enddo

 ! write w' field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'WP 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) wp(i,k,j)
   enddo
  enddo
 enddo
endif



if(vorflag.eq.1)then
 ! write omx field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'OMx 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) omx(i,k,j)
   enddo
  enddo
 enddo

 ! write omy field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'OMy 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) omy(i,k,j)
   enddo
  enddo
 enddo

 ! write omz field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'OMz 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) omz(i,k,j)
   enddo
  enddo
 enddo
endif



if(strflag.eq.1)then
 ! write str x field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'STRAINx 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) strx(i,k,j)
   enddo
  enddo
 enddo

 ! write str y field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'STRAINy 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) stry(i,k,j)
   enddo
  enddo
 enddo

 ! write str z field
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer = 'STRAINz 1 '//str4//' double'//lf
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  do j=y_start,y_end,dny
   do i=x_start,x_end,dnx
    write(66) strz(i,k,j)
   enddo
  enddo
 enddo
endif



 ! write phi field
 if(phiflag.eq.1)then
  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'PHI 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=z_start,z_end,dnz
   do j=y_start,y_end,dny
    do i=x_start,x_end,dnx
     write(66) phi(i,k,j)
    enddo
   enddo
  enddo
 endif

 ! write psi field
 if(psiflag.eq.1)then
  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'PSI 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=z_start,z_end,dnz
   do j=y_start,y_end,dny
    do i=x_start,x_end,dnx
     write(66) psi(i,k,j)
    enddo
   enddo
  enddo
 endif

 ! write theta field
 if(tempflag.eq.1)then
  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'T 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=z_start,z_end,dnz
   do j=y_start,y_end,dny
    do i=x_start,x_end,dnx
     write(66) theta(i,k,j)
    enddo
   enddo
  enddo
 endif

 buffer=lf
 write(66) trim(buffer)

 close(66,status='keep')

return
end
