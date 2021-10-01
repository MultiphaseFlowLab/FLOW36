program main

use commondata 

integer :: start,last,step
integer :: i,k

start=4000000
last =4000000
step =10000


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

allocate(greuu(nz),greuv(nz),greuw(nz),grevv(nz),grevw(nz),greww(nz))
greuu=0.0d0
greuv=0.0d0
greuw=0.0d0
grevv=0.0d0
grevw=0.0d0
greww=0.0d0
count=0

! start averaging from st_count
st_count=0

!call create_plan

do i=start,last,step
  call calc_shear(i)
enddo

!call destroy_plan
!gtau=gtau/dble(count)
greuu=greuu/dble(count)
greuv=greuv/dble(count)
greuw=greuw/dble(count)
grevv=grevv/dble(count)
grevw=grevw/dble(count)
greww=greww/dble(count)

open(546,file='./output/shear.dat',form='formatted',status='new')
write(546,'(7a16))') 'z','re_uu','re_uv','re_uw','re_vv','re_vw','re_ww'
do k=1,nz
  write(546,'(7es16.6))') dcos((dble(k-1)*3.14159265358979)/dble(nz-1)),greuu(k),greuv(k),greuw(k),grevv(k),grevw(k),greww(k)
enddo
close(546,status='keep')


deallocate(greuu,greuv,greuw,grevv,grevw,greww)

return
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_shear(nstep)

use commondata 

integer :: nstep
integer :: i,j,k

double precision, dimension(nx,nz,ny) :: u,v,w
double precision, dimension(nx,nz,ny) :: re_uu,re_uv,re_uw,re_vv,re_vw,re_ww
!double precision, dimension(nx/2+1,nz,ny,2) :: uc,tauc
double precision, dimension(nz) :: um,vm,wm
double precision, dimension(nz) :: re_uum,re_uvm,re_uwm,re_vvm,re_vwm,re_wwm

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

re_uu=0.0d0
re_uv=0.0d0
re_uw=0.0d0
re_vv=0.0d0
re_vw=0.0d0
re_ww=0.0d0

do j=1,ny
 do k=1,nz
  do i=1,nx
    re_uu(i,k,j)=(u(i,k,j)-um(k))*(u(i,k,j)-um(k))
    re_uv(i,k,j)=(u(i,k,j)-um(k))*(v(i,k,j)-vm(k))
    re_uw(i,k,j)=(u(i,k,j)-um(k))*(w(i,k,j)-wm(k))
    re_vv(i,k,j)=(v(i,k,j)-vm(k))*(v(i,k,j)-vm(k))
    re_vw(i,k,j)=(v(i,k,j)-vm(k))*(w(i,k,j)-wm(k))
    re_ww(i,k,j)=(w(i,k,j)-wm(k))*(w(i,k,j)-wm(k))
  enddo
 enddo
enddo


re_uum=0.0d0
re_uvm=0.0d0
re_uwm=0.0d0
re_vvm=0.0d0
re_vwm=0.0d0
re_wwm=0.0d0

do j=1,ny
 do k=1,nz
  do i=1,nx
    re_uum(k)=re_uum(k)+re_uu(i,k,j)
    re_uvm(k)=re_uvm(k)+re_uv(i,k,j)
    re_uwm(k)=re_uwm(k)+re_uw(i,k,j)
    re_vvm(k)=re_vvm(k)+re_vv(i,k,j)
    re_vwm(k)=re_vwm(k)+re_vw(i,k,j)
    re_wwm(k)=re_wwm(k)+re_ww(i,k,j)
  enddo
 enddo
enddo

re_uum=re_uum/dble(nx*ny)
re_uvm=re_uvm/dble(nx*ny)
re_uwm=re_uwm/dble(nx*ny)
re_vvm=re_vvm/dble(nx*ny)
re_vwm=re_vwm/dble(nx*ny)
re_wwm=re_wwm/dble(nx*ny)


if(nstep.ge.st_count)then
 count=count+1
 greuu=greuu+re_uum
 greuv=greuv+re_uvm
 greuw=greuw+re_uwm
 grevv=grevv+re_vvm
 grevw=grevw+re_vwm
 greww=greww+re_wwm
endif


fname='./output/shear_'//time//'.dat'
open(34,file=trim(fname),status='new',form='formatted')
write(34,'(7(a16))') 'z','re_uu','re_uv','re_uw','re_vv','re_vw','re_ww'
do k=1,nz
  write(34,'(7(es16.6))') dcos((dble(k-1)*3.14159265358979)/dble(nz-1)),re_uum(k),re_uvm(k),re_uwm(k),re_vvm(k),re_vwm(k),re_wwm(k)
enddo
close(34,status='keep')


return
end

