subroutine courant_check(u,v,w)

use commondata
use mpi
use sim_par
use par_size
use grid
use phase_field
use surfactant

double precision, dimension(nx,fpz,fpy) :: u,v,w
double precision :: lcomax,gcomax

integer :: i,j,k

#define phi_flag phicompflag
#define psi_flag psicompflag
#define machine machineflag

lcomax=0.0d0
do j=2,fpy
  do k=2,fpz
    do i=2,nx
      lcomax=max(lcomax,dt*(dabs(u(i-1,k,j)/(x(i)-x(i-1)))+ &
         &         dabs(v(i,k,j-1)/(y(fstart(3)+j)-y(fstart(3)+j-1)))+ &
         &         dabs(w(i,k-1,j)/(z(fstart(2)+k)-z(fstart(2)+k-1)))))
    enddo
  enddo
enddo

#if phi_flag == 1
#if machine == 4
if(isnan(phic(1,1,1,1)).eq.1) lcomax=7.0d0
#else
if(isnan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#endif

#if psi_flag == 1
#if machine == 4
if(isnan(psic_fg(1,1,1,1)).eq.1) lcomax=7.0d0
#else
if(isnan(psic_fg(1,1,1,1)).eqv..true.) lcomax=7.0d0
#endif
#endif
#endif

call mpi_allreduce(lcomax,gcomax,1,mpi_double,mpi_max,mpi_comm_world,ierr)

if(rank.eq.0) write(*,'(1x,a,es8.2)') 'check on Courant number : ',gcomax

#if machine == 4
if((gcomax.gt.co).or.(isnan(gcomax).eq.1))then
  if(rank.eq.0) write(*,'(1x,a,es8.2,a,f8.3)') 'Courant number exceeds maximum value : ',gcomax,'>',co
  call exit(0)
endif
#else
if((gcomax.gt.co).or.(isnan(gcomax).eqv..true.))then
  if(rank.eq.0) write(*,'(1x,a,es8.2,a,f8.3)') 'Courant number exceeds maximum value : ',gcomax,'>',co
  call exit(0)
endif
#endif

return
end
