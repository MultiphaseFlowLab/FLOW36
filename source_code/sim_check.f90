subroutine sim_check(i,int_1)

use commondata
use sim_par
use velocity
use phase_field
use surfactant
use dual_grid

integer :: i

double precision :: time,re_bulk,int_phi,int_1,int_psi
double precision :: tau_wm1,tau_wp1,nu_wm1,nu_wp1

#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define cpiflag cpicompflag

time=re*dble(i)*dt

! re_bulk=2.0d0*re/dble(2*nx*ny)*((uc(1,1,1,1))**2+(vc(1,1,1,1))**2+(wc(1,1,1,1)**2))**0.5d0
call get_rebulk(re_bulk)

#if phiflag == 1
call phi_check(int_phi)
#if psiflag == 1
call psi_check(int_psi)
#endif
#endif

#if tempflag == 1
call theta_check(tau_wm1,tau_wp1,nu_wm1,nu_wp1)
#endif

open(66,file=trim(folder)//'/time_check.dat',status='old',position='append')

#if phiflag == 0 && tempflag == 0 && cpiflag == 0
 write(66,'(2(2x,es16.5))') time, re_bulk
#elif phiflag == 0 && tempflag == 0 && cpiflag == 1
 write(66,'(3(2x,es16.5))') time, re_bulk, gradpx
#elif phiflag == 1 && tempflag == 0 && psiflag == 0 && cpiflag == 0
 write(66,'(4(2x,es16.5))') time, re_bulk, int_phi, int_1
#elif phiflag == 1 && tempflag == 0 && psiflag == 0 && cpiflag == 1
 write(66,'(5(2x,es16.5))') time, re_bulk, gradpx, int_phi, int_1
#elif phiflag == 1 && tempflag == 0 && psiflag == 1
  write(66,'(5(2x,es16.5))') time, re_bulk, int_phi, int_1, int_psi
#elif phiflag == 0 && tempflag == 1
 write(66,'(6(2x,es16.5))') time,re_bulk,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#elif phiflag == 1 && tempflag == 1 && psiflag == 0
 write(66,'(8(2x,es16.5))') time, re_bulk, int_phi, int_1,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#elif phiflag == 1 && tempflag == 1 && psiflag == 1
 write(66,'(9(2x,es16.5))') time, re_bulk, int_phi, int_1,int_psi,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#endif

close(66,status='keep')


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_check(int_1)

use commondata
use sim_par
use velocity
use phase_field
use surfactant
use dual_grid

double precision :: time,re_bulk,int_phi,int_1,int_psi
double precision :: tau_wm1,tau_wp1,nu_wm1,nu_wp1

#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag

time=re*0.0d0

! re_bulk=2.0d0*re/dble(2*nx*ny)*((uc(1,1,1,1))**2+(vc(1,1,1,1))**2+(wc(1,1,1,1)**2))**0.5d0
call get_rebulk(re_bulk)

#if phiflag == 1
call phi_check(int_phi)
#if psiflag == 1
call psi_check(int_psi)
#endif
#endif

#if tempflag == 1
call theta_check(tau_wm1,tau_wp1,nu_wm1,nu_wp1)
#endif

open(66,file=trim(folder)//'/time_check.dat',status='new')
open(67,file=trim(folder)//'/backup/time_check.dat',status='new')
open(68,file=trim(folder)//'/backup/time_check_old.dat',status='new')

#if phiflag == 0 && tempflag == 0 && cpiflag == 0
 write(66,'(2(2x,a16))') 't+ ','Re_bulk '
 write(66,'(2(2x,es16.5))') time, re_bulk
 write(67,'(2(2x,a16))') 't+ ','Re_bulk '
 write(67,'(2(2x,es16.5))') time, re_bulk
#elif phiflag == 0 && tempflag == 0 && cpiflag == 1
 write(66,'(3(2x,a16))') 't+ ','Re_bulk ','Gradpx'
 write(66,'(3(2x,es16.5))') time, re_bulk, gradpx
 write(67,'(3(2x,a16))') 't+ ','Re_bulk ','Gradpx'
 write(67,'(3(2x,es16.5))') time, re_bulk, gradpx
#elif phiflag == 1 && tempflag == 0 && psiflag == 0 && cpiflag == 0
 write(66,'(4(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1'
 write(66,'(4(2x,es16.5))') time, re_bulk,int_phi,int_1
 write(67,'(4(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1'
 write(67,'(4(2x,es16.5))') time, re_bulk,int_phi,int_1
#elif phiflag == 1 && tempflag == 0 && psiflag == 0 && cpiflag == 1
 write(66,'(5(2x,a16))') 't+ ','Re_bulk ','Gradpx','phi avg','int phi=+1'
 write(66,'(5(2x,es16.5))') time, re_bulk, gradpx, int_phi,int_1
 write(67,'(5(2x,a16))') 't+ ','Re_bulk ','Gradpx','phi avg','int phi=+1'
 write(67,'(5(2x,es16.5))') time, re_bulk, gradpx, int_phi,int_1
#elif phiflag == 1 && tempflag == 0 && psiflag == 1
 write(66,'(5(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1','psi avg'
 write(66,'(5(2x,es16.5))') time, re_bulk,int_phi,int_1,int_psi
 write(67,'(5(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1','psi avg'
 write(67,'(5(2x,es16.5))') time, re_bulk,int_phi,int_1,int_psi
#elif phiflag == 0 && tempflag == 1
 write(66,'(6(2x,a16))') 't+ ','Re_bulk ','tau wall at z=-1','tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(66,'(6(2x,es16.5))') time,re_bulk,tau_wm1,tau_wp1,nu_wm1,nu_wp1
 write(67,'(6(2x,a16))') 't+ ','Re_bulk ','tau wall at z=-1','tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(67,'(6(2x,es16.5))') time,re_bulk,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#elif phiflag == 1 && tempflag == 1 && psiflag == 0
 write(66,'(8(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1','tau wall at z=-1', &
 &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(66,'(8(2x,es16.5))') time, re_bulk, int_phi, int_1,tau_wm1,tau_wp1,nu_wm1,nu_wp1
 write(67,'(8(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1','tau wall at z=-1', &
 &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(67,'(8(2x,es16.5))') time, re_bulk, int_phi, int_1,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#elif phiflag == 1 && tempflag == 1 && psiflag == 1
  write(66,'(9(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1','psi avg','tau wall at z=-1', &
  &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
  write(66,'(9(2x,es16.5))') time, re_bulk, int_phi, int_1,int_psi,tau_wm1,tau_wp1,nu_wm1,nu_wp1
  write(67,'(9(2x,a16))') 't+ ','Re_bulk ','phi avg','int phi=+1','psi avg','tau wall at z=-1', &
  &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
  write(67,'(9(2x,es16.5))') time, re_bulk, int_phi, int_1,int_psi,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#endif

close(66,status='keep')
close(67,status='keep')
close(68,status='keep')

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine integral_phi(int_1)

use commondata
use grid
use par_size
use phase_field
use mpi

double precision :: int_1, sum_xy(fpz),sglob(nz),column(nz)
double precision :: top(nx,fpz,fpy)

integer :: k

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      if(phi(i,k,j).gt.0.0d0)then
        top(i,k,j)=1.0d0
      else
        top(i,k,j)=0.0d0
      endif
    enddo
  enddo
enddo

sum_xy=sum(sum(top,3),1)

sglob=0.0d0
sglob(fstart(2)+1:fstart(2)+fpz)=sum_xy

call mpi_reduce(sglob,column,nz,mpi_double_precision,mpi_sum,0,flow_comm,ierr)

int_1=0.0d0
! integrate column and save in int_1
do k=1,nz-1
  int_1=int_1+(column(k)+column(k+1))*(z(k)-z(k+1))*0.5d0
enddo

! int_1 *dx*dy
int_1=int_1*xl*yl/dble((nx-1)*(ny-1))

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine theta_check(tau_wm1,tau_wp1,nu_wm1,nu_wp1)

use commondata
use par_size
use velocity
use temperature

double precision, dimension(spx,nz,spy,2) :: tau,nu
double precision :: tau_c(1,nz,1,2),nu_c(1,nz,1,2)
double precision :: tau_wm1,tau_wp1,nu_wm1,nu_wp1

call dz(uc+vc,tau)

tau_c=0.0d0
tau_c(1,:,1,1)=tau(1,:,1,1)

call dctz_bwd(tau_c,tau_c,1,nz,1,0)

tau_wm1=tau_c(1,nz,1,1)/dble(nx*ny)
tau_wp1=tau_c(1,1,1,1)/dble(nx*ny)


call dz(thetac,nu)

nu_c=0.0d0
nu_c(1,:,1,1)=nu(1,:,1,1)

call dctz_bwd(nu_c,nu_c,1,nz,1,0)

nu_wm1=nu_c(1,nz,1,1)/dble(nx*ny)
nu_wp1=nu_c(1,1,1,1)/dble(nx*ny)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phi_check(int_phi)

use commondata
use par_size
use phase_field
use grid

double precision :: int_phi
double precision, dimension(nz) :: dz,phi_mm
double precision, dimension(nz,2) :: phi_mmc

integer :: k


dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do i=2,nz-1
  dz(i)=(z(i-1)-z(i+1))*0.5d0
enddo

call dctz_bwd(phic(1,:,1,1),phi_mmc,1,nz,1,0)

phi_mm=phi_mmc(:,1)/dble(nx*ny)

int_phi=0.0d0
do k=1,nz
  int_phi=int_phi+phi_mm(k)*dz(k)
enddo

! average phi
int_phi=int_phi/2.0d0

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psi_check(int_psi)

use commondata
use par_size
use surfactant
use grid

double precision :: int_psi
double precision, dimension(nz) :: dz,psi_mm
double precision, dimension(nz,2) :: psi_mmc

integer :: k


dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do i=2,nz-1
  dz(i)=(z(i-1)-z(i+1))*0.5d0
enddo

call dctz_bwd(psic(1,:,1,1),psi_mmc,1,nz,1,0)

psi_mm=psi_mmc(:,1)/dble(nx*ny)

int_psi=0.0d0
do k=1,nz
  int_psi=int_psi+psi_mm(k)*dz(k)
enddo

! average psi
int_psi=int_psi/2.0d0

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_rebulk(re_bulk)

use commondata
use par_size
use velocity
use grid
use sim_par

double precision :: re_bulk
double precision, dimension(nz) :: dz,ubulk,vbulk,wbulk
double precision, dimension(nz,2) :: tmp_mmc

integer :: k


dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do i=2,nz-1
  dz(i)=(z(i-1)-z(i+1))*0.5d0
enddo

call dctz_bwd(uc(1,:,1,1),tmp_mmc,1,nz,1,0)
ubulk=tmp_mmc(:,1)/dble(nx*ny)

call dctz_bwd(vc(1,:,1,1),tmp_mmc,1,nz,1,0)
vbulk=tmp_mmc(:,1)/dble(nx*ny)

call dctz_bwd(wc(1,:,1,1),tmp_mmc,1,nz,1,0)
wbulk=tmp_mmc(:,1)/dble(nx*ny)


re_bulk=0.0d0
! sums up Re_b,x+Re_b,y+Re_b,z
do k=1,nz
  re_bulk=re_bulk+(ubulk(k)+vbulk(k)+wbulk(k))*dz(k)
enddo

! average velocity (module)
re_bulk=re_bulk/2.0d0

! bulk Reynolds number (channel height)
re_bulk=re_bulk*re

#if cpiflag == 1
  gradpx=Repow*re/(2*re_bulk)
  print*,'New gradpx',gradpx
#endif


return
end
