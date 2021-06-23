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
    write(*,'(2x,a,i5,a)') 'New drop, ',drop_count,' drops'
    ! single drop part
    s_drop=0
    s_drop(id,kd,jd)=1
    ! flood fill algorithm
    call flood_fill(top,s_drop,id,jd,kd)
    ! remove drops already done from top
    top=top-s_drop
    ! calculate drop volume and center of mass
    call get_center(s_drop,nstep)
    ! new drop calculation
   endif
  enddo
 enddo
enddo

write(*,'(2x,a,i5)') 'Number of drops: ',drop_count
write(*,*)

open(42,file='./output/drop_count.dat',access='append',form='formatted',status='old')
 write(42,'(i16,2x,es16.6,2x,i16)') nstep,dble(nstep)*dt*re,drop_count
close(42,status='keep')

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_center(s_drop,nstep)

use commondata

double precision :: sum_xy(nz)
double precision :: int_1,dx,dy,dz
double precision :: xg,yg,zg
double precision, dimension(nx,nz,ny) :: phix,phiy,phiz

integer :: s_drop(nx,nz,ny)
integer :: nstep
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

call normal(phix,phiy,phiz)

call princ_axis_2d(s_drop,xg,yg,zg,nstep,xboundary,yboundary)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine normal(phix,phiy,phiz)

use commondata
use wavenumbers

double precision, dimension(nx/2+1,nz,ny,2) :: sphix,sphiy,sphiz,tphic
double precision, dimension(nx,nz,ny) :: phix,phiy,phiz,tphi
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

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine princ_axis_2d(s_drop,xg,yg,zg,nstep,xboundary,yboundary)

use commondata

double precision :: xg,yg,zg
double precision :: a11,a12,a22,cff
double precision :: dx,dy,dz,xit,yjt
double precision, dimension(2) :: eigvalue,eigv1,eigv2
double precision :: a(2,2)
double precision :: a1,a2,L,B,D,alpha

integer :: s_drop(nx,nz,ny),nstep,xboundary,yboundary
integer :: i,j,k


a11=0.0d0
a12=0.0d0
a22=0.0d0

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)


! for 3D simulations
do j=ny/2,ny/2
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
   a11=a11+dble(s_drop(i,k,j))*(x(i)-xg)**2*dx*dz
   a12=a12-dble(s_drop(i,k,j))*(x(i)-xg)*(z(k)-zg)*dx*dz
   a22=a22+dble(s_drop(i,k,j))*(z(k)-zg)**2*dx*dz
  enddo
 enddo
enddo


a(1,1)=a11
a(1,2)=a12
a(2,1)=a12
a(2,2)=a22

eigvalue(1)=(-(-a(1,1)-a(2,2))+sqrt((-a(1,1)-a(2,2))**2-4.0d0*(a(1,1)*a(2,2)-a(1,2)*a(2,1))))/2.0d0;
eigvalue(2)=(-(-a(1,1)-a(2,2))-sqrt((-a(1,1)-a(2,2))**2-4.0d0*(a(1,1)*a(2,2)-a(1,2)*a(2,1))))/2.0d0;

a2=((4.0d0*eigvalue(1)/pi)**3*pi/(4.0d0*eigvalue(2)))**(1.0/8.0)
a1=4.0d0*eigvalue(1)/(pi*a2**3)


! first eigenvector
eigv1(1)=1.0d0
eigv1(2)=-a(2,1)/(a(2,2)-eigvalue(1))

!normalize by cff
cff=dsqrt(1+(a(2,1)/(a(2,2)-eigvalue(1)))**2)
eigv1=eigv1/cff


! second eigenvector
eigv2(1)=1.0d0
eigv2(2)=-a(2,1)/(a(2,2)-eigvalue(2))

!normalize by cff
cff=dsqrt(1+(a(2,1)/(a(2,2)-eigvalue(2)))**2)
eigv2=eigv2/cff


! calculate inclination in radians
if(eigvalue(1).ge.eigvalue(2))then
  alpha=atan(-eigv1(2)/eigv1(1))
else
  alpha=atan(-eigv2(2)/eigv2(1))
endif

! output angle in degrees
alpha=alpha*180.0d0/pi

! Taylor deformation
L=max(a1,a2)
B=min(a1,a2)
D=(L-B)/(L+B)

!write(*,*) L,B,D


!! Tryggvason deformation sqrt(Imax/Imin)
!L=max(eigvalue(1),eigvalue(2))
!B=min(eigvalue(1),eigvalue(2))
!D=sqrt(L/B)


open(55,file='./output/deformation.dat',form='formatted',access='append',status='old')
write(55,'(i8,es16.5,2(es16.8))') nstep,dt*re*dble(nstep),D,alpha
close(55,status='keep')


return
end

