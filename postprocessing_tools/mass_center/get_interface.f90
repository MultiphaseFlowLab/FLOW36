subroutine get_interface(nstep)

use commondata

integer :: nstep
integer :: i,j,k,id,jd,kd
integer :: top(nx,nz,ny),s_drop(nx,nz,ny)
integer :: drop_count


do j=1,ny
  do k=1,nz
    do i=1,nx
      if(phi(i,k,j).ge.0.0d0)then
        top(i,k,j)=1
      else
        top(i,k,j)=0
      endif
    enddo
  enddo
enddo

drop_count=0

! flood fill algorithm
do jd=1,ny
 do kd=1,nz
  do id=1,nx
   if(top(id,kd,jd).gt.0)then
    drop_count=drop_count+1
    write(*,'(2x,a,i3,a)') 'New drop, ',drop_count,' drops'
    ! single drop part
    s_drop=0
    s_drop(id,kd,jd)=1
    ! flood fill algorithm
    call flood_fill(top,s_drop,id,jd,kd)
    ! remove drops already done from top
    top=top-s_drop
    ! calculate drop volume and center of mass
    call get_center(s_drop,nstep,drop_count)
    ! new drop calculation
   endif
  enddo
 enddo
enddo

write(*,'(2x,a,i4)') 'Number of drops: ',drop_count
write(*,*)

open(42,file='./output/drop_count.dat',access='append',form='formatted',status='old')
 write(42,'(i16,2x,es16.6,2x,i16)') nstep,dble(nstep)*dt*re,drop_count
close(42,status='keep')

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_center(s_drop,nstep,drop_count)

use commondata

double precision :: sum_xy(nz)!,inertia_t(6)
double precision :: int_1,dx,dy,dz
double precision :: xg,yg,zg,A!,Mx,My,Mz

integer :: s_drop(nx,nz,ny)
integer :: nstep,drop_count
integer :: i,j,k,xboundary,yboundary,it,jt

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)
! dx, dy are uniform in space, dz=dz(z), integral of phi=+1 on V can be reduced to
! a sum of integrals over z
sum_xy=sum(sum(dble(s_drop),3),1)

int_1=0.0d0;
do k=1,nz-1
  int_1=int_1+(sum_xy(k)+sum_xy(k+1))*(z(k)-z(k+1))*0.5d0
enddo
int_1=int_1*dx*dy

xboundary=0
yboundary=0

do j=1,ny
  do k=1,nz
    do i=1,nx
      if(s_drop(i,k,j).eq.1)then
       if((i.eq.nx).or.(i.eq.1)) xboundary=1
       if((j.eq.ny).or.(j.eq.1)) yboundary=1
      endif
    enddo
  enddo
enddo

xg=0.0d0
yg=0.0d0
zg=0.0d0
do j=1,ny
  do k=1,nz-1
    dz=z(k)-z(k+1)
    do i=1,nx
      if(xboundary.eq.1) then
       it=mod(i-floor(real(nx)/4.0)+nx,nx)+1
      else
       it=i
      endif
      if(yboundary.eq.1) then
       jt=mod(j-floor(real(ny)/4.0)+ny,ny)+1
      else
       jt=j
      endif
      xg=xg+dble(s_drop(i,k,j))*dx*dy*dz*x(it)
      yg=yg+dble(s_drop(i,k,j))*dx*dy*dz*y(jt)
      zg=zg+dble(s_drop(i,k,j))*dx*dy*dz*z(k)
    enddo
  enddo
enddo

xg=xg/int_1
yg=yg/int_1
zg=zg/int_1

if(xboundary.eq.1)then
 if(floor(real(nx)/4.0).ge.1)then
  xg=xg+x(floor(real(nx)/4.0))
  if(xg.gt.xl) xg=xg-xl
 endif
endif

if(yboundary.eq.1)then
  if(floor(real(ny)/4.0).ge.1)then
   yg=yg+y(floor(real(ny)/4.0))
   if(yg.gt.yl) yg=yg-yl
  endif
endif


call get_area(s_drop,A)
!call get_area(s_drop,A,xg,yg,zg,Mx,My,Mz,xboundary,yboundary)

!call princ_axis(s_drop,xg,yg,zg,inertia_t,nstep,drop_count,xboundary,yboundary)

! comparison with volume of a sphere
! write(*,'(i8,5(f12.6))') nstep,xg,yg,zg,int_1,4.0d0/3.0d0*pi*0.3d0**3

open(66,file='./output/mass_center.dat',form='formatted',status='old',access='append')
! output in outer units
 ! write(66,'(i8,i6,5(es16.6))') nstep,drop_count,xg,yg,zg,int_1,A
! output in wall units
! write(66,'(i8,i6,8(es16.6))') nstep,drop_count,xg*Re,yg*Re,(zg+1.0d0)*Re,int_1*Re**3,A*Re**2,Mx*Re**3,My*Re**3,Mz*Re**3
 write(66,'(i8,i6,5(es16.6))') nstep,drop_count,xg*Re,yg*Re,(zg+1.0d0)*Re,int_1*Re**3,A*Re**2
close(66,status='keep')

!open(67,file='./output/inertia_tensor.dat',form='formatted',status='old',access='append')
! write(67,'(2(i6),6(es16.6))') nstep,drop_count,inertia_t(1:6)
!close(67,status='keep')

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_area(s_drop,A)
! subroutine get_area(s_drop,A,xg,yg,zg,Mx,My,Mz,xboundary,yboundary)

use commondata

double precision :: A!,norm(3),space(3),Aij,rx,ry,rz,dMx,dMy,dMz,xg,yg,zg
!double precision :: Mx,My,Mz,xit,yjt
double precision :: dx,dy,dz(nz),threshold

integer :: s_drop(nx,nz,ny)
integer :: i,j,k

! logical :: answer



A=0.0d0
threshold=0.5d0
! Mx=0.0d0
! My=0.0d0
! Mz=0.0d0

dz(1)=(z(1)-z(2))*0.5d0
dz(nz)=(z(nz-1)-z(nz))*0.5d0
do k=2,nz-1
 dz(k)=(z(k-1)-z(k+1))*0.5d0
enddo

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)
do j=1,ny
 do k=1,nz
  do i=1,nx
   ! if(s_drop(i,k,j).eq.1) A=A+div(i,k,j)*dx*dy*dz(k)
   if((s_drop(i,k,j).eq.1).and.(phi(i,k,j).lt.dtanh(threshold*ch/(dsqrt(2.0d0)*Ch))))then
    A=A+dx*dy*dz(k)
   endif
  enddo
 enddo
enddo

A=A/(threshold*ch)

! space(1)=xl/dble(nx-1)
! space(3)=yl/dble(ny-1)
! do j=1,ny
!  do k=2,nz-1
!   do i=1,nx
!    norm(1)=phix(i,k,j)
!    norm(2)=phiz(i,k,j)
!    norm(3)=phiy(i,k,j)
!    space(2)=(z(k-1)-z(k+1))*0.5d0
!    call check_interface(answer,i,j,k,s_drop)
!    if(answer.eqv..true.)then
!    ! if(abs(phi(i,k,j)).lt.0.02)then
!     call get_area_ij(norm,space,Aij)
!     A=A+abs(Aij)
!     ! only for torque calculation
!     ! check if drop is crossing a boundary: xboundary, yboundary
!     ! xit=x(i)
!     ! yjt=y(j)
!     ! if(xboundary.eq.1)then
!     !  if(xg.lt.0.5d0*xl)then
!     !    if(xit.ge.0.5d0*xl) xit=xit-xl
!     !  else
!     !    if(xit.lt.0.5d0*xl) xit=xit+xl
!     !  endif
!     ! endif
!     ! if(yboundary.eq.1)then
!     !  if(yg.lt.0.5d0*yl)then
!     !    if(yjt.ge.0.5d0*yl) yjt=yjt-yl
!     !  else
!     !    if(yjt.lt.0.5d0*yl) yjt=yjt+yl
!     !  endif
!     ! endif
!     ! rx=xit-xg
!     ! ry=yjt-yg
!     ! rz=z(k)-zg
!     ! call torque(rx,ry,rz,dmx,dmy,dmz,i,j,k)
!     ! Mx=Mx+dmx*space(1)*space(2)*space(3)
!     ! My=My+dmy*space(1)*space(2)*space(3)
!     ! Mz=Mz+dmz*space(1)*space(2)*space(3)
!    endif
!   enddo
!  enddo
! enddo


! write(*,*) A,4.0d0*pi*0.5**2

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_area_ij(norm,space,area)

double precision, dimension(3) :: norm,space
double precision :: area

double precision, dimension(6,3) :: inter
double precision :: a,b,c,d
double precision :: dx,dy,dz
double precision :: xp,yp,zp,xc,yc,zc
double precision :: c1,c2,c3

integer :: int_count,i,i1,i2

! equation of plane: a*x+b*y+c*z=d
a=norm(1)
b=norm(3)
c=norm(2)

dx=space(1)
dy=space(3)
dz=space(2)

! temp only
xc=dx*0.5d0
yc=dy*0.5d0
zc=dz*0.5d0
! find point with phi=0




! d set so that plane crosses center of cell
d=a*xc+b*yc+c*zc


! find all intersections plane with cell, from 3 to 6
int_count=0
! loop over all 12 edges
yp=0.0d0
zp=0.0d0
xp=(d-b*yp-c*zp)/a
if((xp.ge.0.0d0).and.(xp.le.dx))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

yp=dy
zp=0.0d0
xp=(d-b*yp-c*zp)/a
if((xp.ge.0.0d0).and.(xp.le.dx))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

yp=0.0d0
zp=dz
xp=(d-b*yp-c*zp)/a
if((xp.ge.0.0d0).and.(xp.le.dx))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

yp=dy
zp=dz
xp=(d-b*yp-c*zp)/a
if((xp.ge.0.0d0).and.(xp.le.dx))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif



xp=0.0d0
zp=0.0d0
yp=(d-a*xp-c*zp)/b
if((yp.ge.0.0d0).and.(yp.le.dy))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

xp=dx
zp=0.0d0
yp=(d-a*xp-c*zp)/b
if((yp.ge.0.0d0).and.(yp.le.dy))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

xp=0.0d0
zp=dz
yp=(d-a*xp-c*zp)/b
if((yp.ge.0.0d0).and.(yp.le.dy))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

xp=dx
zp=dz
yp=(d-a*xp-c*zp)/b
if((yp.ge.0.0d0).and.(yp.le.dy))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif



xp=0.0d0
yp=0.0d0
zp=(d-a*xp-b*yp)/c
if((zp.ge.0.0d0).and.(zp.le.dz))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

xp=dx
yp=0.0d0
zp=(d-a*xp-b*yp)/c
if((zp.ge.0.0d0).and.(zp.le.dz))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

xp=0.0d0
yp=dy
zp=(d-a*xp-b*yp)/c
if((zp.ge.0.0d0).and.(zp.le.dz))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

xp=dx
yp=dy
zp=(d-a*xp-b*yp)/c
if((zp.ge.0.0d0).and.(zp.le.dz))then
  ! intersection found
  int_count=int_count+1
  inter(int_count,1)=xp
  inter(int_count,2)=yp
  inter(int_count,3)=zp
endif

if((int_count>6).or.(int_count<3)) write(*,*) 'count intersections (>=3, <=6): ',int_count


area=0.0d0
! calculate area of triangles making the polygon
do i=1,int_count
!write(*,*) i,mod(i,int_count)+1
  i1=i
  i2=mod(i,int_count)+1
  c1=(inter(i2,2)-inter(i1,2))*(zc-inter(i1,3))-(inter(i2,3)-inter(i1,3))*(yc-inter(i1,2))
  c2=(inter(i2,3)-inter(i1,3))*(xc-inter(i1,1))-(inter(i2,1)-inter(i1,1))*(zc-inter(i1,3))
  c3=(inter(i2,1)-inter(i1,1))*(yc-inter(i1,2))-(inter(i2,2)-inter(i1,2))*(xc-inter(i1,1))
  area=area+abs(dsqrt(c1**2+c2**2+c3**2))
enddo

area=area*0.5d0

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine normal

use commondata
use wavenumbers

double precision, dimension(nx/2+1,nz,ny,2) :: sphix,sphiy,sphiz,tphic!,divc
double precision, dimension(nx,nz,ny) :: tphi
double precision :: modphi
integer :: i,j,k

! outward pointing normal
tphi=-phi

call phys_to_spectral(tphi,tphic,1)


! calculate normal to the interface
do j=1,ny
 do k=1,nz
  do i=1,nx/2+1
   sphix(i,k,j,1)=-kx(i)*tphic(i,k,j,2)
   sphix(i,k,j,2)=kx(i)*tphic(i,k,j,1)
   sphiy(i,k,j,1)=-ky(j)*tphic(i,k,j,2)
   sphiy(i,k,j,2)=ky(j)*tphic(i,k,j,1)
  enddo
 enddo
enddo

call dz(tphic,sphiz)


call spectral_to_phys(sphix,phix,1)
call spectral_to_phys(sphiy,phiy,1)
call spectral_to_phys(sphiz,phiz,1)


do j=1,ny
 do k=1,nz-1
  do i=1,nx
   modphi=sqrt(phix(i,k,j)**2+phiy(i,k,j)**2+phiz(i,k,j)**2)
   phix(i,k,j)=phix(i,k,j)/modphi
   phiy(i,k,j)=phiy(i,k,j)/modphi
   phiz(i,k,j)=phiz(i,k,j)/modphi
 enddo
enddo
enddo

! call phys_to_spectral(phix,sphix,0)
! call phys_to_spectral(phiy,sphiy,0)
! call phys_to_spectral(phiz,sphiz,0)
!
! ! calculate divergence
! call dz(sphiz,divc)
! do j=1,ny
!  do k=1,nz
!   do i=1,nx/2+1
!    divc(i,k,j,1)=divc(i,k,j,1)-kx(i)*sphix(i,k,j,2)-ky(j)*sphiy(i,k,j,2)
!    divc(i,k,j,2)=divc(i,k,j,2)+kx(i)*sphix(i,k,j,1)+ky(j)*sphiy(i,k,j,1)
!   enddo
!  enddo
! enddo
!
! call spectral_to_phys(divc,div,0)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_interface(answer,i,j,k,s_drop)

use commondata

logical :: answer,boundary

integer :: s_drop(nx,nz,ny)

integer :: i,j,k
integer :: ii,jj

boundary=.false.


! could be done with elseif
if(s_drop(i,k,j).eq.1)then
 ! check orthogonal direction
 !xdir
 ii=mod(i,nx)+1
 if(s_drop(ii,k,j).eq.0)then
  boundary=.true.
 endif
 ii=mod(i+nx-2,nx)+1
 if(s_drop(ii,k,j).eq.0)then
  boundary=.true.
 endif
! ydir
 jj=mod(j,ny)+1
 if(s_drop(i,k,jj).eq.0)then
  boundary=.true.
 endif
 jj=mod(j+ny-2,ny)+1
 if(s_drop(i,k,jj).eq.0)then
  boundary=.true.
 endif
! zdir
 if(s_drop(i,k+1,j).eq.0)then
  boundary=.true.
 endif
 if(s_drop(i,k-1,j).eq.0)then
  boundary=.true.
 endif
! check on edges (all 3x3x3 cube)
 ! i+1,j+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ii=mod(i,nx)+1
 jj=mod(j,ny)+1
 if(s_drop(ii,k-1,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k+1,jj).eq.0)then
  boundary=.true.
 endif
 ! i+1,j-1
 jj=mod(j+ny-2,ny)+1
 if(s_drop(ii,k-1,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k+1,jj).eq.0)then
  boundary=.true.
 endif
 ! i+1,j
 if(s_drop(ii,k-1,j).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k,j).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k+1,j).eq.0)then
  boundary=.true.
 endif
 ! i,j+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 jj=mod(j,ny)+1
 if(s_drop(i,k-1,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(i,k,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(i,k+1,jj).eq.0)then
  boundary=.true.
 endif
 ! i,j-1
 jj=mod(j+ny-2,ny)+1
 if(s_drop(i,k-1,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(i,k,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(i,k+1,jj).eq.0)then
  boundary=.true.
 endif
 ! i,j
 if(s_drop(i,k-1,j).eq.0)then
  boundary=.true.
 endif
 ! if(s_drop(i,k,j).eq.0)then
 !  boundary=.true.
 ! endif
 if(s_drop(i,k+1,j).eq.0)then
  boundary=.true.
 endif
 ! i-1,j+1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ii=mod(i+nx-2,nx)+1
 jj=mod(j,ny)+1
 if(s_drop(ii,k-1,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k+1,jj).eq.0)then
  boundary=.true.
 endif
 ! i-1,j-1
 jj=mod(j+ny-2,ny)+1
 if(s_drop(ii,k-1,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k,jj).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k+1,jj).eq.0)then
  boundary=.true.
 endif
 ! i-1,j
 if(s_drop(ii,k-1,j).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k,j).eq.0)then
  boundary=.true.
 endif
 if(s_drop(ii,k+1,j).eq.0)then
  boundary=.true.
 endif

endif


! if at least one neighbour is equal to zero and the node value is 1 it is an interface node
if((s_drop(i,k,j).eq.1).and.(boundary.eqv..true.))then
 answer=.true.
else
 answer=.false.
endif

! if
! if((s_drop(i,k,j).eq.1).and.())then
!  answer=.true.
! else
!  answer=.false.
! endif


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine princ_axis(s_drop,xg,yg,zg,inertia_t,nstep,count,xboundary,yboundary)

use commondata

double precision :: xg,yg,zg
double precision :: a11,a12,a13,a22,a23,a33,cff,inertia_t(6)
double precision :: dx,dy,dz,xit,yjt
double precision, dimension(3) :: eigv1,eigv2,eigv3,eigvalue
double precision :: a(3,3)

integer :: s_drop(nx,nz,ny),nstep,count,xboundary,yboundary
integer :: i,j,k


a11=0.0d0
a12=0.0d0
a13=0.0d0
a22=0.0d0
a23=0.0d0
a33=0.0d0

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)


do j=1,ny
 do k=2,nz-1
  dz=(z(k-1)-z(k+1))/2.0d0
  do i=1,nx
   xit=x(i)
   yjt=y(j)
   if(xboundary.eq.1)then
    if(xg.lt.0.5d0*xl)then
      if(xit.ge.0.5d0*xl) xit=xit-xl
    else
      if(xit.lt.0.5d0*xl) xit=xit+xl
    endif
   endif
   if(yboundary.eq.1)then
    if(yg.lt.0.5d0*yl)then
      if(yjt.ge.0.5d0*yl) yjt=yjt-yl
    else
      if(yjt.lt.0.5d0*yl) yjt=yjt+yl
    endif
   endif
   a11=a11+s_drop(i,k,j)*((yjt-yg)**2+(z(k)-zg)**2)*dx*dy*dz
   a12=a12-s_drop(i,k,j)*((xit-xg)*(yjt-yg))*dx*dy*dz
   a13=a13-s_drop(i,k,j)*((xit-xg)*(z(k)-zg))*dx*dy*dz
   a22=a22+s_drop(i,k,j)*((xit-xg)**2+(z(k)-zg)**2)*dx*dy*dz
   a23=a23-s_drop(i,k,j)*((yjt-yg)*(z(k)-zg))*dx*dy*dz
   a33=a33+s_drop(i,k,j)*((xit-xg)**2+(yjt-yg)**2)*dx*dy*dz
  enddo
 enddo
enddo

a(1,1)=a11
a(1,2)=a12
a(1,3)=a13
a(2,1)=a(1,2)
a(2,2)=a22
a(2,3)=a23
a(3,1)=a(1,3)
a(3,2)=a(2,3)
a(3,3)=a33

! transform matrix A into upper triangular matrix
! eliminate 1st column
cff=a(2,1)/a(1,1)
a(2,1)=a(2,1)-cff*a(1,1)
a(2,2)=a(2,2)-cff*a(1,2)
a(2,3)=a(2,3)-cff*a(1,3)
cff=a(3,1)/a(1,1)
a(3,1)=a(3,1)-cff*a(1,1)
a(3,2)=a(3,2)-cff*a(1,2)
a(3,3)=a(3,3)-cff*a(1,3)

! eliminate 2nd column
cff=a(3,2)/a(2,2)
a(3,2)=a(3,2)-cff*a(2,2)
a(3,3)=a(3,3)-cff*a(2,3)

! do i=1,3
!  write(*,'(3(es20.8))') a(i,:)
! enddo


eigvalue(1)=a(1,1)
eigvalue(2)=a(2,2)
eigvalue(3)=a(3,3)

! first eigenvector
eigv1(1)=1.0d0
eigv1(2)=0.0d0
eigv1(3)=0.0d0

! second eigenvector + normalization
eigv2(1)=-a(1,2)/(eigvalue(1)-eigvalue(2))
eigv2(2)=1.0d0
eigv2(3)=0.0d0
cff=sqrt((eigv2(1))**2+(eigv2(2))**2+(eigv2(3))**2)
eigv2=eigv2/cff

! second eigenvector + normalization
eigv3(1)=a(2,3)*a(1,2)/(eigvalue(2)-eigvalue(3))-a(1,3)/(eigvalue(1)-eigvalue(3))
eigv3(2)=-a(2,3)/(eigvalue(2)-eigvalue(3))
eigv3(3)=1.0d0
cff=sqrt((eigv3(1))**2+(eigv3(2))**2+(eigv3(3))**2)
eigv3=eigv3/cff


! do i=1,3
!  write(*,'(3(es8.1,2x))')  eigv1(i),eigv2(i),eigv3(i)
! enddo

!open(68,file='./output/eigen_problem.dat',form='formatted',access='append',status='old')
! write(68,'(2(i6),4(es16.6))') nstep,count,eigvalue(1),eigv1
! write(68,'(2(i6),4(es16.6))') nstep,count,eigvalue(2),eigv2
! write(68,'(2(i6),4(es16.6))') nstep,count,eigvalue(3),eigv3
!close(68,status='keep')


! to be removed
inertia_t(1)=a11
inertia_t(2)=a12
inertia_t(3)=a13
inertia_t(4)=a22
inertia_t(5)=a23
inertia_t(6)=a33



return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine torque(rx,ry,rz,dmx,dmy,dmz,i,j,k)

use sterm

double precision :: rx,ry,rz,dmx,dmy,dmz

integer :: i,j,k

dmx=ry*s3(i,k,j)-rz*s2(i,k,j)
dmy=rz*s1(i,k,j)-rx*s3(i,k,j)
dmz=rx*s2(i,k,j)-ry*s1(i,k,j)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
