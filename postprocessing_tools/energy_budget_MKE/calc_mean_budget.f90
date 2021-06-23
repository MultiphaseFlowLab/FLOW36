subroutine calc_budget_mean(step)

use commondata
use fields
use wavenumber

integer :: step
integer :: i,j,k

double precision, dimension(nz) :: volume
double precision, dimension(nz) :: um,vm,wm
double precision, dimension(nx,nz,ny) :: rho,eta
double precision, dimension(nxfg,nzfg,nyfg) :: sigma
double precision, dimension(nxfg/2+1,nzfg,nyfg,2) :: phicfg
double precision, dimension(nz) :: press_work,visc_diff,visc_diss,prod,turb_diff,interf,gravity,balance
double precision, allocatable, dimension(:) :: a1,a1c
double precision, allocatable, dimension(:,:,:) :: a11,a12,a13,a22,a23,a33
double precision, allocatable, dimension(:,:,:,:) :: a11c,a12c,a13c,a22c,a23c,a33c,a1cc,a2cc,a3cc
double precision, allocatable, dimension(:,:,:) :: phix,phiy,phiz


double precision :: dx,dy,V_int,V_d,V_c,int_thick
double precision :: int_pwork_c,int_vdiff_c,int_vdiss_c,int_prod_c,int_tdiff_c,int_gravity_c
double precision :: int_pwork_d,int_vdiff_d,int_vdiss_d,int_prod_d,int_tdiff_d,int_gravity_d,int_interf
character(len=40) :: namefile

! write(*,*) 'Computing MKE budgets.......'

write(*,'(1x,a,i8,a,i8)') 'Step ',step,' of ',nend

! write(*,*) 'Starting to computed MKE....'
! write(*,*) 'Reading fields..............'

call read_fields(step)

! thickness of interface
int_thick=dtanh(3.0d0*ch/(dsqrt(2.0d0)*Ch))

! volume of interface, droplets and carrier
V_int=0.0d0
V_d=0.0d0
V_c=0.0d0
! ! integral quantities
int_pwork_c=0.0d0
int_vdiff_c=0.0d0
int_vdiss_c=0.0d0
int_prod_c=0.0d0
int_tdiff_c=0.0d0
int_gravity_c=0.0d0
int_pwork_d=0.0d0
int_vdiff_d=0.0d0
int_vdiss_d=0.0d0
int_prod_d=0.0d0
int_tdiff_d=0.0d0
int_gravity_d=0.0d0
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pressure diffusion
! write(*,*) 'Pressure work...............'
press_work=0.0d0
press_work=-gradpx*um-gradpy*vm

do j=1,ny
 do k=1,nz
  do i=1,nx
   if(phi(i,k,j).ge.0.0d0)then
    V_d=V_d+volume(k)
    int_pwork_d=int_pwork_d+volume(k)*press_work(k)
   else
    V_c=V_c+volume(k)
    int_pwork_c=int_pwork_c+volume(k)*press_work(k)
   endif
  enddo
 enddo
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(*,*) 'Viscous diffusion...........'
allocate(a1(nz))
allocate(a1c(nz))

a1=um**2.0d0

call dctz_fwd_1d(a1,a1c,nz,0)

call dz_1d(a1c,a1)
call dz_1d(a1,a1c)

call dctz_bwd_1d(a1c,a1,nz,0)

visc_diff=0.5d0/re*a1

deallocate(a1,a1c)

do j=1,ny
 do k=1,nz
  do i=1,nx
   if(phi(i,k,j).ge.0.0d0)then
    int_vdiff_d=int_vdiff_d+volume(k)*visc_diff(k)
   else
    int_vdiff_c=int_vdiff_c+volume(k)*visc_diff(k)
   endif
  enddo
 enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(*,*) 'Viscous dissipation.........'
visc_diss=0.0d0

allocate(a1(nz))
allocate(a1c(nz))

call dctz_fwd_1d(um,a1c,nz,0)

call dz_1d(a1c,a1)

call dctz_bwd_1d(a1,a1c,nz,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   visc_diss(k)=visc_diss(k)+eta(i,k,j)*(a1c(k)*a1c(k))
  enddo
 enddo
enddo

deallocate(a1,a1c)

visc_diss=-1.0d0/re*visc_diss/dble(nx*ny)

do j=1,ny
 do k=1,nz
  do i=1,nx
   if(phi(i,k,j).ge.0.0d0)then
    int_vdiss_d=int_vdiss_d+volume(k)*visc_diss(k)
   else
    int_vdiss_c=int_vdiss_c+volume(k)*visc_diss(k)
   endif
  enddo
 enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(*,*) 'TKE production..............'
prod=0.0d0

allocate(a11(nx,nz,ny))
allocate(a33(nx,nz,ny))

allocate(a1(nz))
allocate(a1c(nz))

do j=1,ny
 do i=1,nx
  a11(i,:,j)=u(i,:,j)-um
  a33(i,:,j)=w(i,:,j)-wm
 enddo
enddo

call dctz_fwd_1d(um,a1,nz,0)

call dz_1d(a1,a1c)

call dctz_bwd_1d(a1c,a1,nz,0)

do j=1,ny
 do k=1,nz
  do i=1,nx
   prod(k)=prod(k)+rho(i,k,j)*(a1(k)*a11(i,k,j)*a33(i,k,j))
   if(phi(i,k,j).ge.0.0d0)then
    int_prod_d=int_prod_d+volume(k)*rho(i,k,j)*(a1(k)*a11(i,k,j)*a33(i,k,j))
   else
    int_prod_c=int_prod_c+volume(k)*rho(i,k,j)*(a1(k)*a11(i,k,j)*a33(i,k,j))
   endif
  enddo
 enddo
enddo

prod=prod/dble(nx*ny)

deallocate(a1,a1c)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(*,*) 'Turbulent transport.........'
turb_diff=0.0d0

allocate(a11c(nx/2+1,nz,ny,2))
allocate(a33c(nx/2+1,nz,ny,2))

do j=1,ny
 do k=1,nz
  do i=1,nx
   a11(i,k,j)=um(k)*a11(i,k,j)*a33(i,k,j)
  enddo
 enddo
enddo

deallocate(a33)

call phys_to_spectral(a11,a11c,0)

call dz(a11c,a33c)

deallocate(a11c)

call spectral_to_phys(a33c,a11,0)

deallocate(a33c)

do j=1,ny
 do k=1,nz
  do i=1,nx
   turb_diff(k)=turb_diff(k)+rho(i,k,j)*a11(i,k,j)
   if(phi(i,k,j).ge.0.0d0)then
    int_tdiff_d=int_tdiff_d+volume(k)*rho(i,k,j)*a11(i,k,j)
   else
    int_tdiff_c=int_tdiff_c+volume(k)*rho(i,k,j)*a11(i,k,j)
   endif
  enddo
 enddo
enddo

deallocate(a11)

turb_diff=-turb_diff/dble(nx*ny)
int_tdiff_c=-int_tdiff_c
int_tdiff_d=-int_tdiff_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write(*,*) 'Interfacial term............'
interf=0.0d0
! leave to zero if no interface
if(phi_flag.eq.1)then
  ! port phi to fine grid
  call phys_to_spectral(phi,phic,0)
  call coarse2fine(phic,phicfg)

  ! assemble Korteweg tensor in fine grid and multiply by EOS surface tension

  allocate(a11(nxfg,nzfg,nyfg))
  allocate(a12(nxfg,nzfg,nyfg))
  allocate(a13(nxfg,nzfg,nyfg))
  allocate(a22(nxfg,nzfg,nyfg))
  allocate(a23(nxfg,nzfg,nyfg))
  allocate(a33(nxfg,nzfg,nyfg))
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

  allocate(phix(nxfg,nzfg,nyfg))
  allocate(phiy(nxfg,nzfg,nyfg))
  allocate(phiz(nxfg,nzfg,nyfg))

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

  deallocate(phix,phiy,phiz)

  ! take divergence of tensor
  call phys_to_spectral_fg(a11,a11c,0)
  call phys_to_spectral_fg(a12,a12c,0)
  call phys_to_spectral_fg(a13,a13c,0)
  call phys_to_spectral_fg(a22,a22c,0)
  call phys_to_spectral_fg(a23,a23c,0)
  call phys_to_spectral_fg(a33,a33c,0)

  deallocate(a11,a12,a13,a22,a23,a33)

  allocate(a1cc(nxfg/2+1,nzfg,nyfg,2))
  allocate(a2cc(nxfg/2+1,nzfg,nyfg,2))
  allocate(a3cc(nxfg/2+1,nzfg,nyfg,2))

  call dz_fg(a13c,a1cc)
  call dz_fg(a23c,a2cc)
  call dz_fg(a33c,a3cc)
  do j=1,nyfg
   do i=1,nxfg/2+1
     a1cc(i,:,j,1)=a1cc(i,:,j,1)-kxfg(i)*a11c(i,:,j,2)-kyfg(j)*a12c(i,:,j,2)
     a1cc(i,:,j,2)=a1cc(i,:,j,2)+kxfg(i)*a11c(i,:,j,1)+kyfg(j)*a12c(i,:,j,1)
     a2cc(i,:,j,1)=a2cc(i,:,j,1)-kxfg(i)*a12c(i,:,j,2)-kyfg(j)*a22c(i,:,j,2)
     a2cc(i,:,j,2)=a2cc(i,:,j,2)+kxfg(i)*a12c(i,:,j,1)+kyfg(j)*a22c(i,:,j,1)
     a3cc(i,:,j,1)=a3cc(i,:,j,1)-kxfg(i)*a13c(i,:,j,2)-kyfg(j)*a23c(i,:,j,2)
     a3cc(i,:,j,2)=a3cc(i,:,j,2)+kxfg(i)*a13c(i,:,j,1)+kyfg(j)*a23c(i,:,j,1)
   enddo
  enddo

  ! port to coarse grid
  deallocate(a11c,a12c,a13c,a22c,a23c,a33c)
  allocate(a11c(nx/2+1,nz,ny,2))
  allocate(a22c(nx/2+1,nz,ny,2))
  allocate(a33c(nx/2+1,nz,ny,2))

  call fine2coarse(a1cc,a11c)
  call fine2coarse(a2cc,a22c)
  call fine2coarse(a3cc,a33c)

  deallocate(a1cc,a2cc,a3cc)

  allocate(a11(nx,nz,ny))
  allocate(a22(nx,nz,ny))
  allocate(a33(nx,nz,ny))

  call spectral_to_phys(a11c,a11,0)
  call spectral_to_phys(a22c,a22,0)
  call spectral_to_phys(a33c,a33,0)

  deallocate(a11c,a22c,a33c)

  do j=1,ny
   do k=1,nz
    do i=1,nx
     interf(k)=interf(k)+um(k)*a11(i,k,j)+vm(k)*a22(i,k,j)+wm(k)*a33(i,k,j)
     if(dabs(phi(i,k,j)).le.int_thick)then
      V_int=V_int+volume(k)
      int_interf=int_interf+volume(k)*(um(k)*a11(i,k,j)+vm(k)*a22(i,k,j)+wm(k)*a33(i,k,j))
     endif
    enddo
   enddo
  enddo

  deallocate(a11,a22,a33)

  interf=3.0d0/dsqrt(8.0d0)*Ch/We*interf/dble(nx*ny)
  int_interf=3.0d0/dsqrt(8.0d0)*Ch/We*int_interf
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gravity term
gravity=0.0d0
if(phi_flag.eq.1)then
 if(b_type.eq.0)then ! no buoyancy and gravity
  gravity=0.0d0
 elseif(b_type.eq.1)then ! buoyancy and gravity (rho*g)
  do j=1,ny
   do k=1,nz
    do i=1,nx
     gravity(k)=gravity(k)+rho(i,k,j)*(um(k)*grav(1)+vm(k)*grav(3)+wm(k)*grav(2))/Fr**2.0d0
     if(phi(i,k,j).ge.0.0d0)then
      int_gravity_d=int_gravity_d+volume(k)*rho(i,k,j)*(um(k)*grav(1)+vm(k)*grav(3)+wm(k)*grav(2))/Fr**2.0d0
     else
      int_gravity_c=int_gravity_c+volume(k)*rho(i,k,j)*(um(k)*grav(1)+vm(k)*grav(3)+wm(k)*grav(2))/Fr**2.0d0
     endif
    enddo
   enddo
  enddo
 elseif(b_type.eq.2)then ! only buoyancy (Delta rho*g)
  do j=1,ny
   do k=1,nz
    do i=1,nx
     gravity(k)=gravity(k)+(rho(i,k,j)-1.0d0)*(um(k)*grav(1)+vm(k)*grav(3)+wm(k)*grav(2))/Fr**2.0d0
     if(phi(i,k,j).ge.0.0d0)then
      int_gravity_d=int_gravity_d+volume(k)*(rho(i,k,j)-1.0d0)*(um(k)*grav(1)+vm(k)*grav(3)+wm(k)*grav(2))/Fr**2.0d0
     else
      int_gravity_c=int_gravity_c+volume(k)*(rho(i,k,j)-1.0d0)*(um(k)*grav(1)+vm(k)*grav(3)+wm(k)*grav(2))/Fr**2.0d0
     endif
    enddo
   enddo
  enddo
 endif
endif

gravity=gravity/dble(nx*ny)




! normalize energy budget in w.u.
press_work=press_work/re
visc_diff=visc_diff/re
visc_diss=visc_diss/re
prod=prod/re
turb_diff=turb_diff/re
interf=interf/re
gravity=gravity/re


int_pwork_c=int_pwork_c/re
int_vdiff_c=int_vdiff_c/re
int_vdiss_c=int_vdiss_c/re
int_prod_c=int_prod_c/re
int_tdiff_c=int_tdiff_c/re
int_gravity_c=int_gravity_c/re

int_pwork_d=int_pwork_d/re
int_vdiff_d=int_vdiff_d/re
int_vdiss_d=int_vdiss_d/re
int_prod_d=int_prod_d/re
int_tdiff_d=int_tdiff_d/re
int_gravity_d=int_gravity_d/re

int_interf=int_interf/re


balance=press_work+visc_diff+visc_diss+prod+turb_diff+interf+gravity

! write MKE budget to file
write(namefile,'(a,i8.8,a)') './output/mkebudget_',step,'.dat'

open(54,file=namefile,status='new',form='formatted')
write(54,'(9(a16))') 'z^+','press work','visc diff','visc diss','prod','turb diff','interf','gravity','balance'
! to get from 0 to 2*Re
do k=nz,1,-1
 write(54,'(9(e16.5))') z(k),press_work(k),visc_diff(k),visc_diss(k),prod(k),turb_diff(k),interf(k),gravity(k),balance(k)
enddo
close(54,status='keep')

! write integral budgets to file
open(123,file='./output/integral_mean.dat',status='old',form='formatted',access='append')
write(123,'(i10,14(es15.6))') step,dble(step)*re*dt,int_pwork_c,int_vdiff_c,int_vdiss_c,int_prod_c,int_tdiff_c,int_gravity_c, &
&                                 int_pwork_d,int_vdiff_d,int_vdiss_d,int_prod_d,int_tdiff_d,int_gravity_d,int_interf
close(123,status='keep')


return
end
