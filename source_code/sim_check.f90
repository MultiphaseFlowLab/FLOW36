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

time=re*dble(i)*dt

re_bulk=2.0d0*re/dble(2*nx*ny)*((uc(1,1,1,1))**2+(vc(1,1,1,1))**2+(wc(1,1,1,1)**2))**0.5d0

#if phiflag == 1
int_phi=phic(1,1,1,1)/dble(2*nx*ny)
#if psiflag == 1
int_psi=psic_fg(1,1,1,1)/dble(2*npsix*npsiy)
#endif
#endif

#if tempflag == 1
call theta_check(tau_wm1,tau_wp1,nu_wm1,nu_wp1)
#endif

open(66,file=trim(folder)//'/time_check.dat',status='old',position='append')

#if phiflag == 0 && tempflag == 0
 write(66,'(2(2x,es16.5))') time, re_bulk
#elif phiflag == 1 && tempflag == 0 && psiflag == 0
 write(66,'(4(2x,es16.5))') time, re_bulk, int_phi, int_1
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

re_bulk=2.0d0*re/dble(2*nx*ny)*((uc(1,1,1,1))**2+(vc(1,1,1,1))**2+(wc(1,1,1,1)**2))**0.5d0

#if phiflag == 1
int_phi=phic(1,1,1,1)/dble(2*nx*ny)
#if psiflag == 1
int_psi=psic_fg(1,1,1,1)/dble(2*npsix*npsiy)
#endif
#endif

#if tempflag == 1
call theta_check(tau_wm1,tau_wp1,nu_wm1,nu_wp1)
#endif

open(66,file=trim(folder)//'/time_check.dat',status='new')
open(67,file=trim(folder)//'/backup/time_check.dat',status='new')
open(68,file=trim(folder)//'/backup/time_check_old.dat',status='new')

#if phiflag == 0 && tempflag == 0
 write(66,'(2(2x,a16))') 't+ ','Re_bulk '
 write(66,'(2(2x,es16.5))') time, re_bulk
 write(67,'(2(2x,a16))') 't+ ','Re_bulk '
 write(67,'(2(2x,es16.5))') time, re_bulk
#elif phiflag == 1 && tempflag == 0 && psiflag == 0
 write(66,'(4(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1'
 write(66,'(4(2x,es16.5))') time, re_bulk,int_phi,int_1
 write(67,'(4(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1'
 write(67,'(4(2x,es16.5))') time, re_bulk,int_phi,int_1
#elif phiflag == 1 && tempflag == 0 && psiflag == 1
 write(66,'(5(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1','int psi'
 write(66,'(5(2x,es16.5))') time, re_bulk,int_phi,int_1,int_psi
 write(67,'(5(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1','int psi'
 write(67,'(5(2x,es16.5))') time, re_bulk,int_phi,int_1,int_psi
#elif phiflag == 0 && tempflag == 1
 write(66,'(6(2x,a16))') 't+ ','Re_bulk ','tau wall at z=-1','tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(66,'(6(2x,es16.5))') time,re_bulk,tau_wm1,tau_wp1,nu_wm1,nu_wp1
 write(67,'(6(2x,a16))') 't+ ','Re_bulk ','tau wall at z=-1','tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(67,'(6(2x,es16.5))') time,re_bulk,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#elif phiflag == 1 && tempflag == 1 && psiflag == 0
 write(66,'(8(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1','tau wall at z=-1', &
 &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(66,'(8(2x,es16.5))') time, re_bulk, int_phi, int_1,tau_wm1,tau_wp1,nu_wm1,nu_wp1
 write(67,'(8(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1','tau wall at z=-1', &
 &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
 write(67,'(8(2x,es16.5))') time, re_bulk, int_phi, int_1,tau_wm1,tau_wp1,nu_wm1,nu_wp1
#elif phiflag == 1 && tempflag == 1 && psiflag == 1
  write(66,'(9(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1','int psi','tau wall at z=-1', &
  &                           'tau wall at z=+1','Nu at z=-1','Nu at z=+1'
  write(66,'(9(2x,es16.5))') time, re_bulk, int_phi, int_1,int_psi,tau_wm1,tau_wp1,nu_wm1,nu_wp1
  write(67,'(9(2x,a16))') 't+ ','Re_bulk ','int phi','int phi=+1','int psi','tau wall at z=-1', &
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

call mpi_reduce(sglob,column,nz,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

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
