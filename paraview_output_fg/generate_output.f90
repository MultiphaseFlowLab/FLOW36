subroutine generate_output(nstep)

use commondata
use wavenumber

integer :: nstep,nfields,numx,numy,numz
integer :: i,k,j

double precision :: meanu,meanv,meanw

double precision, dimension(nxf/2+1,nzf,nyf,2) :: tmpc

character(len=40) :: namefile
character(len=80) :: buffer
character(len=10) :: lf
character(len=8) :: str1,str2,str3
character(len=16) :: str4

! end of line character
lf=achar(10)

! fields included
nfields=uflag+vflag+wflag+phiflag+psiflag+tempflag+3*upflag

! calculate velocity fluctuations
if(upflag.eq.1)then
  do k=1,nzf
    meanu=0.0d0
    meanv=0.0d0
    meanw=0.0d0
    do j=1,nyf
      do i=1,nxf
        meanu=meanu+u(i,k,j)
        meanv=meanv+v(i,k,j)
        meanw=meanw+w(i,k,j)
      enddo
    enddo
    meanu=meanu/dble(nxf*nyf)
    meanv=meanv/dble(nxf*nyf)
    meanw=meanw/dble(nxf*nyf)
    up(:,k,:)=u(:,k,:)-meanu
    vp(:,k,:)=v(:,k,:)-meanv
    wp(:,k,:)=w(:,k,:)-meanw
  enddo
endif

if(vorflag.eq.1)then
! calculate vorticity, u,v,w already available in spectral space
  ! om x
  call dz(-vc,tmpc)
  do j=1,nyf
    do i=1,nxf/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)-ky(j)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)+ky(j)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys_fg(tmpc,omx,0)

  ! om y
  call dz(uc,tmpc)
  do j=1,nyf
    do i=1,nxf/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)+kx(i)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)-kx(i)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys_fg(tmpc,omy,0)

  ! om z
  do j=1,nyf
    do i=1,nxf/2+1
      tmpc(i,:,j,1)=-kx(i)*vc(i,:,j,2)+ky(j)*uc(i,:,j,2)
      tmpc(i,:,j,2)=+kx(i)*vc(i,:,j,1)-ky(j)*uc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys_fg(tmpc,omz,0)
endif


if(strflag.eq.1) then
! calculate strain rate, u,v,w already available in spectral space
  ! str x
  call dz(vc,tmpc)
  do j=1,nyf
    do i=1,nxf/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)-ky(j)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)+ky(j)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,strx,0)
  strx=strx*0.5d0

  ! str y
  call dz(uc,tmpc)
  do j=1,nyf
    do i=1,nxf/2+1
      tmpc(i,:,j,1)=tmpc(i,:,j,1)-kx(i)*wc(i,:,j,2)
      tmpc(i,:,j,2)=tmpc(i,:,j,2)+kx(i)*wc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,stry,0)
  stry=stry*0.5d0

  ! str z
  do j=1,nyf
    do i=1,nxf/2+1
      tmpc(i,:,j,1)=-ky(j)*uc(i,:,j,2)-kx(i)*vc(i,:,j,2)
      tmpc(i,:,j,2)=+ky(j)*uc(i,:,j,1)+kx(i)*vc(i,:,j,1)
    enddo
  enddo
  call spectral_to_phys(tmpc,strz,0)
  strz=strz*0.5d0
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
  write(66) xfg(i)
 enddo
 buffer='Y_COORDINATES '//str2//'  double'//lf ;
 write(66) trim(buffer)
 do j=y_start,y_end,dny
  write(66) yfg(j)
 enddo
 buffer='Z_COORDINATES '//str3//'  double'//lf ;
 write(66) trim(buffer)
 do k=z_start,z_end,dnz
  write(66) zfg(k)
 enddo

 ! write content (data format)
 write(str4(1:16),'(i16)') numx*numy*numz
 buffer='POINT_DATA '//str4//lf
 write(66) trim(buffer)
 write(str1(1:8),'(i8)') nfields
 buffer='FIELD FieldData '//str1//lf
 write(66) trim(buffer)

 ! write u field
if(uflag.eq.1)then
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
endif

if(vflag.eq.1)then
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
endif

if(wflag.eq.1)then
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

 if(upflag.eq.1)then
  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'U_p 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=z_start,z_end,dnz
   do j=y_start,y_end,dny
    do i=x_start,x_end,dnx
     write(66) up(i,k,j)
    enddo
   enddo
  enddo

  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'V_p 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=z_start,z_end,dnz
   do j=y_start,y_end,dny
    do i=x_start,x_end,dnx
     write(66) vp(i,k,j)
    enddo
   enddo
  enddo

  write(str4(1:16),'(i16)') numx*numy*numz
  buffer = 'W_p 1 '//str4//' double'//lf
  write(66) trim(buffer)
  do k=z_start,z_end,dnz
   do j=y_start,y_end,dny
    do i=x_start,x_end,dnx
     write(66) wp(i,k,j)
    enddo
   enddo
  enddo

 endif


 buffer=lf
 write(66) trim(buffer)

 close(66,status='keep')

return
end
