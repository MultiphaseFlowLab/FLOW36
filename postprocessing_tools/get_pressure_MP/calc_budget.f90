subroutine calc_budget(step)

use commondata
use fields
use wavenumber

integer :: step
integer :: i,j,k

double precision, dimension(nz) :: um,vm,wm
double precision, dimension(nx,nz,ny) :: rho,eta
double precision, dimension(nxfg,nzfg,nyfg) :: sigma
double precision, dimension(nxfg/2+1,nzfg,nyfg,2) :: phicfg
double precision, dimension(nz) :: press_diff,visc_diff,visc_diss,prod,turb_diff,interf,balance,volume
double precision, allocatable, dimension(:,:,:) :: a11,a12,a13,a22,a23,a33,a1,a2,a3,phix,phiy,phiz
double precision, allocatable, dimension(:,:,:,:) :: a11c,a12c,a13c,a22c,a23c,a33c,a1c,a2c,a3c
double precision :: dx,dy,V_int,V_d,V_c,int_thick
double precision :: int_pdiff_c,int_vdiff_c,int_vdiss_c,int_prod_c,int_tdiff_c
double precision :: int_pdiff_d,int_vdiff_d,int_vdiss_d,int_prod_d,int_tdiff_d,int_interf
character(len=40) :: namefile

! thickness of interface
int_thick=dtanh(3.0d0*ch/(dsqrt(2.0d0)*Ch))

! volume of interface, droplets and carrier
V_int=0.0d0
V_d=0.0d0
V_c=0.0d0
! integral quantities
int_pdiff_c=0.0d0
int_vdiff_c=0.0d0
int_vdiss_c=0.0d0
int_prod_c=0.0d0
int_tdiff_c=0.0d0
int_pdiff_d=0.0d0
int_vdiff_d=0.0d0
int_vdiss_d=0.0d0
int_prod_d=0.0d0
int_tdiff_d=0.0d0
int_interf=0.0d0


dx=x(2)-x(1)
dy=y(2)-y(1)
volume(1)=dx*dy*(z(1)-z(2))*0.5d0
volume(nz)=dx*dy*(z(nz-1)-z(nz))*0.5d0
do k=2,nz-1
 volume(k)=dx*dy*(z(k-1)-z(k+1))*0.5d0
enddo


if(phi_flag.eq.1)then
  rho=1.0d0+(rhor-1.0d0)/2.0d0*(phi+1.0d0)
  eta=1.0d0+(visr-1.0d0)/2.0d0*(phi+1.0d0)
else
  rho=1.0d0
  eta=1.0d0
endif

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


! u,v,w,press: total fields, mean in um,vm,wm,pm

! remove mean pressure from pressure (average on x and y)
! take space average in x and y


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pressure diffusion
press_diff=0.0d0

allocate(a11(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a33(nx,nz,ny))
allocate(a11c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))
allocate(a12c(nx/2+1,nz,ny,2))

do j=1,ny
 do k=1,nz
  do i=1,nx
   a11(i,k,j)=(u(i,k,j)-um(k))*(press(i,k,j)-pm(k))
   a22(i,k,j)=(v(i,k,j)-vm(k))*(press(i,k,j)-pm(k))
   a33(i,k,j)=(w(i,k,j)-wm(k))*(press(i,k,j)-pm(k))
  enddo
 enddo
enddo

call phys_to_spectral(a11,a11c,0)
call phys_to_spectral(a22,a22c,0)
call phys_to_spectral(a33,a33c,0)

call dz(a33c,a12c)
do j=1,ny
 do i=1,nx/2+1
  a12c(i,:,j,1)=a12c(i,:,j,1)-kx(i)*a11c(i,:,j,2)-ky(j)*a22c(i,:,j,2)
  a12c(i,:,j,2)=a12c(i,:,j,2)+kx(i)*a11c(i,:,j,1)+ky(j)*a22c(i,:,j,1)
 enddo
enddo

call spectral_to_phys(a12c,a11,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   press_diff(k)=press_diff(k)+a11(i,k,j)
   ! compute volume of interface, drop and carrier
   ! compute integral quantities
   if(phi(i,k,j).ge.0)then
    V_d=V_d+volume(k)
    int_pdiff_d=int_pdiff_d+a11(i,k,j)*volume(k)
   else
    V_c=V_c+volume(k)
    int_pdiff_c=int_pdiff_c+a11(i,k,j)*volume(k)
   endif
   if(dabs(phi(i,k,j)).le.int_thick)then
    V_int=V_int+volume(k)
   endif
  enddo
 enddo
enddo
press_diff=-press_diff/dble(nx*ny)
int_pdiff_c=-int_pdiff_c/V_c
! avoid NaN if single phase
if(phi_flag.eq.1) int_pdiff_d=-int_pdiff_d/V_d

deallocate(a11,a22,a33)
deallocate(a11c,a22c,a33c,a12c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! viscous diffusion
visc_diff=0.0d0

allocate(a1(nx,nz,ny))
allocate(a2(nx,nz,ny))
allocate(a3(nx,nz,ny))
allocate(a1c(nx/2+1,nz,ny,2))
allocate(a2c(nx/2+1,nz,ny,2))
allocate(a3c(nx/2+1,nz,ny,2))

do j=1,ny
 do k=1,nz
  do i=1,nx
   a1(i,k,j)=u(i,k,j)-um(k)
   a2(i,k,j)=v(i,k,j)-vm(k)
   a3(i,k,j)=w(i,k,j)-wm(k)
  enddo
 enddo
enddo

! velocity fluctuations
call phys_to_spectral(a1,a1c,0)
call phys_to_spectral(a2,a2c,0)
call phys_to_spectral(a3,a3c,0)

allocate(a11c(nx/2+1,nz,ny,2))
allocate(a12c(nx/2+1,nz,ny,2))
allocate(a13c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a23c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))
allocate(a11(nx,nz,ny))
allocate(a12(nx,nz,ny))
allocate(a13(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a23(nx,nz,ny))
allocate(a33(nx,nz,ny))

call dz(a1c,a13c)
call dz(a2c,a23c)
call dz(a3c,a33c)

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-kx(i)*a1c(i,:,j,2)
  a11c(i,:,j,2)=+kx(i)*a1c(i,:,j,1)
  a22c(i,:,j,1)=-ky(j)*a2c(i,:,j,2)
  a22c(i,:,j,2)=+ky(j)*a2c(i,:,j,1)
  a12c(i,:,j,1)=-kx(i)*a2c(i,:,j,2)-ky(j)*a1c(i,:,j,2)
  a12c(i,:,j,2)=+kx(i)*a2c(i,:,j,1)+ky(j)*a1c(i,:,j,1)
  a13c(i,:,j,1)=a13c(i,:,j,1)-kx(i)*a3c(i,:,j,2)
  a13c(i,:,j,2)=a13c(i,:,j,2)+kx(i)*a3c(i,:,j,1)
  a23c(i,:,j,1)=a23c(i,:,j,1)-ky(j)*a3c(i,:,j,2)
  a23c(i,:,j,2)=a23c(i,:,j,2)+ky(j)*a3c(i,:,j,1)
 enddo
enddo
a12c=0.5d0*a12c
a13c=0.5d0*a13c
a23c=0.5d0*a23c

call spectral_to_phys(a11c,a11,0)
call spectral_to_phys(a12c,a12,0)
call spectral_to_phys(a13c,a13,0)
call spectral_to_phys(a22c,a22,0)
call spectral_to_phys(a23c,a23,0)
call spectral_to_phys(a33c,a33,0)

! calculate S'*u' (physical space)
do j=1,ny
 do k=1,nz
  do i=1,nx
   a11(i,k,j)=a11(i,k,j)*a1(i,k,j)+a12(i,k,j)*a2(i,k,j)+a13(i,k,j)*a3(i,k,j)
   a22(i,k,j)=a12(i,k,j)*a1(i,k,j)+a22(i,k,j)*a2(i,k,j)+a23(i,k,j)*a3(i,k,j)
   a33(i,k,j)=a13(i,k,j)*a1(i,k,j)+a23(i,k,j)*a2(i,k,j)+a33(i,k,j)*a3(i,k,j)
  enddo
 enddo
enddo
! S'*u' (spectral space)
call phys_to_spectral(a11,a11c,0)
call phys_to_spectral(a22,a22c,0)
call phys_to_spectral(a33,a33c,0)

! div(S'*u')
call dz(a33c,a23c)
do j=1,ny
 do i=1,nx/2+1
  a33c(i,:,j,1)=-kx(i)*a11c(i,:,j,2)-ky(j)*a22c(i,:,j,2)+a23c(i,:,j,1)
  a33c(i,:,j,2)=+kx(i)*a11c(i,:,j,1)+ky(j)*a22c(i,:,j,1)+a23c(i,:,j,2)
 enddo
enddo

call spectral_to_phys(a33c,a12,0)

! (nabla eta)*(S'*u')
call phys_to_spectral(eta,a12c,0)

do j=1,ny
 do i=1,nx/2+1
  a13c(i,:,j,1)=-kx(i)*a12c(i,:,j,2)
  a13c(i,:,j,2)=+kx(i)*a12c(i,:,j,1)
  a23c(i,:,j,1)=-ky(j)*a12c(i,:,j,2)
  a23c(i,:,j,2)=+ky(j)*a12c(i,:,j,1)
 enddo
enddo
call spectral_to_phys(a13c,a13,0)
call spectral_to_phys(a23c,a23,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   a12(i,k,j)=eta(i,k,j)*a12(i,k,j)+a13(i,k,j)*a11(i,k,j)+a23(i,k,j)*a22(i,k,j)
  enddo
 enddo
enddo

call dz(a12c,a13c)
call spectral_to_phys(a13c,a13,0)
do j=1,ny
 do k=1,nz
  do i=1,nx
   a12(i,k,j)=a12(i,k,j)+a13(i,k,j)*a33(i,k,j)
  enddo
 enddo
enddo

do j=1,ny
 do k=1,nz
  do i=1,nx
   visc_diff(k)=visc_diff(k)+a12(i,k,j)
   ! compute integral quantities
   if(phi(i,k,j).ge.0)then
    int_vdiff_d=int_vdiff_d+a12(i,k,j)*volume(k)
   else
    int_vdiff_c=int_vdiff_c+a12(i,k,j)*volume(k)
   endif
  enddo
 enddo
enddo

visc_diff=2.0d0/re*visc_diff/dble(nx*ny)
int_vdiff_c=2.0d0/re*int_vdiff_c/V_c
! avoid NaN if single phase
if(phi_flag.eq.1) int_vdiff_d=2.0d0/re*int_vdiff_d/V_d

deallocate(a11,a12,a13,a22,a23,a33)
deallocate(a11c,a12c,a13c,a22c,a23c,a33c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! viscous dissipation
visc_diss=0.0d0

! a1,a2,a3 include u' in physical and spectral space (a1c,a2c,a3c)

allocate(a11c(nx/2+1,nz,ny,2))
allocate(a12c(nx/2+1,nz,ny,2))
allocate(a13c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a23c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))
allocate(a11(nx,nz,ny))
allocate(a12(nx,nz,ny))
allocate(a13(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a23(nx,nz,ny))
allocate(a33(nx,nz,ny))

call dz(a1c,a13c)
call dz(a2c,a23c)
call dz(a3c,a33c)

do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-kx(i)*a1c(i,:,j,2)
  a11c(i,:,j,2)=+kx(i)*a1c(i,:,j,1)
  a22c(i,:,j,1)=-ky(j)*a2c(i,:,j,2)
  a22c(i,:,j,2)=+ky(j)*a2c(i,:,j,1)
  a12c(i,:,j,1)=-kx(i)*a2c(i,:,j,2)-ky(j)*a1c(i,:,j,2)
  a12c(i,:,j,2)=+kx(i)*a2c(i,:,j,1)+ky(j)*a1c(i,:,j,1)
  a13c(i,:,j,1)=a13c(i,:,j,1)-kx(i)*a3c(i,:,j,2)
  a13c(i,:,j,2)=a13c(i,:,j,2)+kx(i)*a3c(i,:,j,1)
  a23c(i,:,j,1)=a23c(i,:,j,1)-ky(j)*a3c(i,:,j,2)
  a23c(i,:,j,2)=a23c(i,:,j,2)+ky(j)*a3c(i,:,j,1)
 enddo
enddo
a12c=0.5d0*a12c
a13c=0.5d0*a13c
a23c=0.5d0*a23c

call spectral_to_phys(a11c,a11,0)
call spectral_to_phys(a12c,a12,0)
call spectral_to_phys(a13c,a13,0)
call spectral_to_phys(a22c,a22,0)
call spectral_to_phys(a23c,a23,0)
call spectral_to_phys(a33c,a33,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   visc_diss(k)=visc_diss(k)+eta(i,k,j)*(a11(i,k,j)*a11(i,k,j)+a22(i,k,j)*a22(i,k,j)+a33(i,k,j)*a33(i,k,j)+ &
 &                    2.0d0*(a12(i,k,j)*a12(i,k,j)+a13(i,k,j)*a13(i,k,j)+a23(i,k,j)*a23(i,k,j)))
   ! compute integral quantities
   if(phi(i,k,j).ge.0)then
    int_vdiss_d=int_vdiss_d+eta(i,k,j)*(a11(i,k,j)*a11(i,k,j)+a22(i,k,j)*a22(i,k,j)+a33(i,k,j)*a33(i,k,j)+ &
  &                    2.0d0*(a12(i,k,j)*a12(i,k,j)+a13(i,k,j)*a13(i,k,j)+a23(i,k,j)*a23(i,k,j)))*volume(k)
   else
    int_vdiss_c=int_vdiss_c+eta(i,k,j)*(a11(i,k,j)*a11(i,k,j)+a22(i,k,j)*a22(i,k,j)+a33(i,k,j)*a33(i,k,j)+ &
  &                    2.0d0*(a12(i,k,j)*a12(i,k,j)+a13(i,k,j)*a13(i,k,j)+a23(i,k,j)*a23(i,k,j)))*volume(k)
   endif
  enddo
 enddo
enddo

visc_diss=-2.0d0/re*visc_diss/dble(nx*ny)
int_vdiss_c=-2.0d0/re*int_vdiss_c/V_c
! avoid NaN if single phase
if(phi_flag.eq.1) int_vdiss_d=-2.0d0/re*int_vdiss_d/V_d

deallocate(a11,a12,a13,a22,a23,a33)
deallocate(a11c,a12c,a13c,a22c,a23c,a33c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! production by mean flow
prod=0.0d0

allocate(a11c(nx/2+1,nz,ny,2))
allocate(a12c(nx/2+1,nz,ny,2))
allocate(a13c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a23c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))
allocate(a11(nx,nz,ny))
allocate(a12(nx,nz,ny))
allocate(a13(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a23(nx,nz,ny))
allocate(a33(nx,nz,ny))

call dz(uc-a1c,a13c)
call dz(vc-a2c,a23c)
call dz(wc-a3c,a33c)

! calculate <S> (mean tensor) from u-u'
do j=1,ny
 do i=1,nx/2+1
  a11c(i,:,j,1)=-kx(i)*(uc(i,:,j,2)-a1c(i,:,j,2))
  a11c(i,:,j,2)=+kx(i)*(uc(i,:,j,1)-a1c(i,:,j,1))
  a22c(i,:,j,1)=-ky(j)*(vc(i,:,j,2)-a2c(i,:,j,2))
  a22c(i,:,j,2)=+ky(j)*(vc(i,:,j,1)-a2c(i,:,j,1))
  a12c(i,:,j,1)=-kx(i)*(vc(i,:,j,2)-a2c(i,:,j,2))-ky(j)*(uc(i,:,j,2)-a1c(i,:,j,2))
  a12c(i,:,j,2)=+kx(i)*(vc(i,:,j,1)-a2c(i,:,j,1))+ky(j)*(uc(i,:,j,1)-a1c(i,:,j,1))
  a13c(i,:,j,1)=a13c(i,:,j,1)-kx(i)*(wc(i,:,j,2)-a3c(i,:,j,2))
  a13c(i,:,j,2)=a13c(i,:,j,2)+kx(i)*(wc(i,:,j,1)-a3c(i,:,j,1))
  a23c(i,:,j,1)=a23c(i,:,j,1)-ky(j)*(wc(i,:,j,2)-a3c(i,:,j,2))
  a23c(i,:,j,2)=a23c(i,:,j,2)+ky(j)*(wc(i,:,j,1)-a3c(i,:,j,1))
 enddo
enddo
a12c=0.5d0*a12c
a13c=0.5d0*a13c
a23c=0.5d0*a23c

call spectral_to_phys(a11c,a11,0)
call spectral_to_phys(a12c,a12,0)
call spectral_to_phys(a13c,a13,0)
call spectral_to_phys(a22c,a22,0)
call spectral_to_phys(a23c,a23,0)
call spectral_to_phys(a33c,a33,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   prod(k)=prod(k)+rho(i,k,j)*(a1(i,k,j)*a1(i,k,j)*a11(i,k,j)+a2(i,k,j)*a2(i,k,j)*a22(i,k,j)+ &
 &     a3(i,k,j)*a3(i,k,j)*a33(i,k,j)+2.0d0*(a1(i,k,j)*a2(i,k,j)*a12(i,k,j)+ &
 &     a1(i,k,j)*a3(i,k,j)*a13(i,k,j)+a2(i,k,j)*a3(i,k,j)*a23(i,k,j)))
   ! compute integral quantities
   if(phi(i,k,j).ge.0)then
    int_prod_d=int_prod_d+rho(i,k,j)*(a1(i,k,j)*a1(i,k,j)*a11(i,k,j)+a2(i,k,j)*a2(i,k,j)*a22(i,k,j)+ &
  &     a3(i,k,j)*a3(i,k,j)*a33(i,k,j)+2.0d0*(a1(i,k,j)*a2(i,k,j)*a12(i,k,j)+ &
  &     a1(i,k,j)*a3(i,k,j)*a13(i,k,j)+a2(i,k,j)*a3(i,k,j)*a23(i,k,j)))*volume(k)
   else
    int_prod_c=int_prod_c+rho(i,k,j)*(a1(i,k,j)*a1(i,k,j)*a11(i,k,j)+a2(i,k,j)*a2(i,k,j)*a22(i,k,j)+ &
  &     a3(i,k,j)*a3(i,k,j)*a33(i,k,j)+2.0d0*(a1(i,k,j)*a2(i,k,j)*a12(i,k,j)+ &
  &     a1(i,k,j)*a3(i,k,j)*a13(i,k,j)+a2(i,k,j)*a3(i,k,j)*a23(i,k,j)))*volume(k)
   endif
  enddo
 enddo
enddo

prod=-prod/dble(nx*ny)
int_prod_c=-int_prod_c/V_c
! avoid NaN if single phase
if(phi_flag.eq.1) int_prod_d=-int_prod_d/V_d

deallocate(a1c,a2c,a3c,a11c,a12c,a13c,a22c,a23c,a33c)
deallocate(a11,a12,a13,a22,a23,a33)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! turbulent diffusion
turb_diff=0.0d0

allocate(a11(nx,nz,ny))
allocate(a22(nx,nz,ny))
allocate(a33(nx,nz,ny))
allocate(a1c(nx/2+1,nz,ny,2))
allocate(a11c(nx/2+1,nz,ny,2))
allocate(a22c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))

do j=1,ny
 do k=1,nz
  do i=1,nx
   a11(i,k,j)=(a1(i,k,j)**2.0d0+a2(i,k,j)**2.0d0+a3(i,k,j)**2.0d0)*a1(i,k,j)
   a22(i,k,j)=(a1(i,k,j)**2.0d0+a2(i,k,j)**2.0d0+a3(i,k,j)**2.0d0)*a2(i,k,j)
   a33(i,k,j)=(a1(i,k,j)**2.0d0+a2(i,k,j)**2.0d0+a3(i,k,j)**2.0d0)*a3(i,k,j)
  enddo
 enddo
enddo

call phys_to_spectral(a11,a11c,0)
call phys_to_spectral(a22,a22c,0)
call phys_to_spectral(a33,a33c,0)

call dz(a33c,a1c)
do j=1,ny
 do i=1,nx/2+1
  a1c(i,:,j,1)=a1c(i,:,j,1)-kx(i)*a11c(i,:,j,2)-ky(j)*a22c(i,:,j,2)
  a1c(i,:,j,2)=a1c(i,:,j,2)+kx(i)*a11c(i,:,j,1)+ky(j)*a22c(i,:,j,1)
 enddo
enddo

call spectral_to_phys(a1c,a33,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   turb_diff(k)=turb_diff(k)+rho(i,k,j)*a33(i,k,j)
   ! compute integral quantities
   if(phi(i,k,j).ge.0)then
    int_tdiff_d=int_tdiff_d+rho(i,k,j)*a33(i,k,j)*volume(k)
   else
    int_tdiff_c=int_tdiff_c+rho(i,k,j)*a33(i,k,j)*volume(k)
   endif
  enddo
 enddo
enddo

turb_diff=-0.5d0*turb_diff/dble(nx*ny)
int_tdiff_c=-0.5d0*int_tdiff_c/V_c
! avoid NaN if single phase
if(phi_flag.eq.1) int_tdiff_d=-0.5d0*int_tdiff_d/V_d

deallocate(a11,a22,a33)
deallocate(a1c,a11c,a22c,a33c)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interface contribution
interf=0.0d0
! leave to zero if no interface
if(phi_flag.eq.1)then
  ! port phi to fine grid
  call coarse2fine(phic,phicfg)

  ! assemble Korteweg tensor in fine grid and multiply by EOS surface tension
  allocate(phix(nxfg,nzfg,nyfg))
  allocate(phiy(nxfg,nzfg,nyfg))
  allocate(phiz(nxfg,nzfg,nyfg))
  allocate(a11(nxfg,nzfg,nyfg))
  allocate(a12(nxfg,nzfg,nyfg))
  allocate(a13(nxfg,nzfg,nyfg))
  allocate(a22(nxfg,nzfg,nyfg))
  allocate(a23(nxfg,nzfg,nyfg))
  allocate(a33(nxfg,nzfg,nyfg))
  allocate(a1c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a2c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a3c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a11c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a12c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a13c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a22c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a23c(nxfg/2+1,nzfg,nyfg,2))
  allocate(a33c(nxfg/2+1,nzfg,nyfg,2))

  call dz_fg(phicfg,a33c)
  do j=1,nyfg
   do i=1,nxfg/2+1
    a11c(i,:,j,1)=-kxfg(i)*phicfg(i,:,j,2)
    a11c(i,:,j,2)=+kxfg(i)*phicfg(i,:,j,1)
    a22c(i,:,j,1)=-kyfg(j)*phicfg(i,:,j,2)
    a22c(i,:,j,2)=+kyfg(j)*phicfg(i,:,j,1)
   enddo
  enddo

  call spectral_to_phys_fg(a11c,phix,0)
  call spectral_to_phys_fg(a22c,phiy,0)
  call spectral_to_phys_fg(a33c,phiz,0)


  ! assemble sigma
  if(psi_flag.eq.1)then
   sigma=max((1.0d0+betas*log(1.0d0-psi)),0.5d0)
  else
   sigma=1.0d0
  endif

  do j=1,nyfg
   do k=1,nzfg
    do i=1,nxfg
     a11(i,k,j)=(phiy(i,k,j)**2.0d0+phiz(i,k,j)**2.0d0)*sigma(i,k,j)
     a12(i,k,j)=(-phix(i,k,j)*phiy(i,k,j))*sigma(i,k,j)
     a13(i,k,j)=(-phix(i,k,j)*phiz(i,k,j))*sigma(i,k,j)
     a22(i,k,j)=(phix(i,k,j)**2.0d0+phiz(i,k,j)**2.0d0)*sigma(i,k,j)
     a23(i,k,j)=(-phiy(i,k,j)*phiz(i,k,j))*sigma(i,k,j)
     a33(i,k,j)=(phix(i,k,j)**2.0d0+phiy(i,k,j)**2.0d0)*sigma(i,k,j)
    enddo
   enddo
  enddo


  ! take divergence of tensor
  call phys_to_spectral_fg(a11,a11c,0)
  call phys_to_spectral_fg(a12,a12c,0)
  call phys_to_spectral_fg(a13,a13c,0)
  call phys_to_spectral_fg(a22,a22c,0)
  call phys_to_spectral_fg(a23,a23c,0)
  call phys_to_spectral_fg(a33,a33c,0)

  call dz_fg(a13c,a1c)
  call dz_fg(a23c,a2c)
  call dz_fg(a33c,a3c)
  do j=1,nyfg
   do i=1,nxfg/2+1
     a1c(i,:,j,1)=a1c(i,:,j,1)-kxfg(i)*a11c(i,:,j,2)-kyfg(j)*a12c(i,:,j,2)
     a1c(i,:,j,2)=a1c(i,:,j,2)+kxfg(i)*a11c(i,:,j,1)+kyfg(j)*a12c(i,:,j,1)
     a2c(i,:,j,1)=a2c(i,:,j,1)-kxfg(i)*a12c(i,:,j,2)-kyfg(j)*a22c(i,:,j,2)
     a2c(i,:,j,2)=a2c(i,:,j,2)+kxfg(i)*a12c(i,:,j,1)+kyfg(j)*a22c(i,:,j,1)
     a3c(i,:,j,1)=a3c(i,:,j,1)-kxfg(i)*a13c(i,:,j,2)-kyfg(j)*a23c(i,:,j,2)
     a3c(i,:,j,2)=a3c(i,:,j,2)+kxfg(i)*a13c(i,:,j,1)+kyfg(j)*a23c(i,:,j,1)
   enddo
  enddo

  ! port to coarse grid
  deallocate(a11c,a22c,a33c)
  allocate(a11c(nx/2+1,nz,ny,2))
  allocate(a22c(nx/2+1,nz,ny,2))
  allocate(a33c(nx/2+1,nz,ny,2))

  call fine2coarse(a1c,a11c)
  call fine2coarse(a2c,a22c)
  call fine2coarse(a3c,a33c)

  call spectral_to_phys(a11c,a11,0)
  call spectral_to_phys(a22c,a22,0)
  call spectral_to_phys(a33c,a33,0)

  ! multiply by u' (scalar product)
  do j=1,ny
   do k=1,nz
    do i=1,nx
     interf(k)=interf(k)+a1(i,k,j)*a11(i,k,j)+a2(i,k,j)*a22(i,k,j)+a3(i,k,j)*a33(i,k,j)
     if(dabs(phi(i,k,j)).le.int_thick)then
      int_interf=int_interf+a1(i,k,j)*a11(i,k,j)+a2(i,k,j)*a22(i,k,j)+a3(i,k,j)*a33(i,k,j)
     endif
    enddo
   enddo
  enddo

  interf=3.0d0/dsqrt(8.0d0)*Ch/We*interf/dble(nx*ny)
  int_interf=3.0d0/dsqrt(8.0d0)*Ch/We*int_interf/V_int

  deallocate(a1,a2,a3)  ! u',v',w' in physical space
  deallocate(a1c,a2c,a3c,a11c,a12c,a13c,a22c,a23c,a33c)
  deallocate(phix,phiy,phiz,a11,a12,a13,a22,a23,a33)
endif

! normalize energy budget in w.u.
press_diff=press_diff/re
visc_diff=visc_diff/re
visc_diss=visc_diss/re
prod=prod/re
turb_diff=turb_diff/re
interf=interf/re

int_pdiff_c=int_pdiff_c/re
int_vdiff_c=int_vdiff_c/re
int_vdiss_c=int_vdiss_c/re
int_prod_c=int_prod_c/re
int_tdiff_c=int_tdiff_c/re
int_pdiff_d=int_pdiff_d/re
int_vdiff_d=int_vdiff_d/re
int_vdiss_d=int_vdiss_d/re
int_prod_d=int_prod_d/re
int_tdiff_d=int_tdiff_d/re
int_interf=int_interf/re

balance=press_diff+visc_diff+visc_diss+prod+turb_diff+interf

! write TKE budget to file
write(namefile,'(a,i8.8,a)') './output/budget_',step,'.dat'

open(54,file=namefile,status='new',form='formatted')

write(54,'(8(a16))') 'z^+','press diff','visc diff','visc diss','prod','turb diff','interf','balance'
! to get from 0 to 2*Re
do k=nz,1,-1
 write(54,'(8(e16.5))') z(k),press_diff(k),visc_diff(k),visc_diss(k),prod(k),turb_diff(k),interf(k),balance(k)
enddo
close(54,status='keep')

! write integral budgets to file
open(123,file='./output/integral.dat',status='old',form='formatted',access='append')
write(123,'(i10,12(es15.6))') step,dble(step)*re*dt,int_pdiff_c,int_vdiff_c,int_vdiss_c,int_prod_c,int_tdiff_c, &
 &                                 int_pdiff_d,int_vdiff_d,int_vdiss_d,int_prod_d,int_tdiff_d,int_interf
close(123,status='keep')


return
end
