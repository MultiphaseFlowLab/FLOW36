program main

use commondata 

integer :: start,last,step
integer :: i,k

start=0
last=55000
step=5000


open(unit=66,file='../sc_compiled/input.f90',form='formatted',status='old',action='read')
read(66,*)
read(66,*)
read(66,*)
read(66,*)
read(66,*)
read(66,'(i5)') nx
read(66,'(i5)') ny
read(66,'(i5)') nz
read(66,*)
read(66,'(f16.6)') Re
close(66,status='keep')

allocate(gtau(nz))
allocate(gre(nz))
gtau=0.0d0
gre=0.0d0
count=0

! start averaging from st_count
st_count=30000

call create_plan

do i=start,last,step
  call calc_shear(i)
enddo

call destroy_plan

gtau=gtau/dble(count)
gre=gre/dble(count)

open(546,file='./output/shear.dat',form='formatted',status='new')
write(546,'(4(a16))') 'z','viscous','re stress','total'
do k=1,nz
  write(546,'(4(es16.6))') dcos((dble(k-1)*3.14159265358979)/dble(nz-1)),gtau(k),gre(k),gtau(k)+gre(k)
enddo
close(546,status='keep')


deallocate(gtau,gre)

return
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_shear(nstep)

use commondata 

integer :: nstep
integer :: i,j,k

double precision, dimension(nx,nz,ny) :: u,v,w
double precision, dimension(nx,nz,ny) :: re_str,tau
double precision, dimension(nx/2+1,nz,ny,2) :: uc,tauc
double precision, dimension(nz) :: um,vm,wm
double precision, dimension(nz) :: taum,restrm

character(len=8) :: time
character(len=160) :: namefile,fname

write(time,'(i8.8)') nstep

namefile='../results/u_'//time//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian') 
read(666) u
close(666,status='keep')

namefile='../results/v_'//time//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) v
close(666,status='keep')

namefile='../results/w_'//time//'.dat'
open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
read(666) w
close(666,status='keep')

write(*,*) nstep, maxval(u),maxval(v),maxval(w)

um=0.0d0
vm=0.0d0
wm=0.0d0

do j=1,ny
 do k=1,nz
  do i=1,nx
   um(k)=um(k)+u(i,k,j)
   vm(k)=vm(k)+v(i,k,j)
   wm(k)=wm(k)+w(i,k,j)
  enddo
 enddo
enddo

um=um/dble(nx*ny)
vm=vm/dble(nx*ny)
wm=wm/dble(nx*ny)


do j=1,ny
 do k=1,nz
  do i=1,nx
    re_str(i,k,j)=(u(i,k,j)-um(k))*(w(i,k,j)-wm(k))
  enddo
 enddo
enddo


call phys_to_spectral(u,uc,0)

call dz(uc,tauc)

call spectral_to_phys(tauc,tau,0)


taum=0.0d0
restrm=0.0d0
do j=1,ny
 do k=1,nz
  do i=1,nx
    taum(k)=taum(k)+tau(i,k,j)
    restrm(k)=restrm(k)+re_str(i,k,j)
  enddo
 enddo
enddo
taum=taum/dble(nx*ny)
! normalization of viscous stress
taum=-taum/re
restrm=restrm/dble(nx*ny)


if(nstep.ge.st_count)then
 count=count+1
 gtau=gtau+taum
 gre=gre+restrm
endif


fname='./output/shear_'//time//'.dat'
open(34,file=trim(fname),status='new',form='formatted')
write(34,'(4(a16))') 'z','viscous','re stress','total'
do k=1,nz
  write(34,'(4(es16.6))') dcos((dble(k-1)*3.14159265358979)/dble(nz-1)),taum(k),restrm(k),taum(k)+restrm(k)
enddo
close(34,status='keep')


return
end

