subroutine courant_check(u,v,w)

use commondata
use mpi
use sim_par
use par_size
use grid
use phase_field
use surfactant
!only for PGI compiler

double precision, dimension(nx,fpz,fpy) :: u,v,w
double precision :: lcomax,gcomax,dx,dy,dz(nz)

integer :: i,j,k

#define phi_flag phicompflag
#define psi_flag psicompflag
#define machine machineflag

!$acc kernels
dx=xl/real(nx-1)
dy=yl/real(ny-1)
dz(1)=z(1)-z(2)
dz(nz)=z(nz-1)-z(nz)
do k=2,nz-1
  dz(k)=(z(k-1)-z(k+1))/2.0d0
enddo
lcomax=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      lcomax=max(lcomax,dt*(dabs(u(i,k,j)/dx)+ &
                            dabs(v(i,k,j)/dy)+ &
                            dabs(w(i,k,j)/dz(fstart(2)+k))))
    enddo
  enddo
enddo
!$acc end kernels

#if machine == 4
if(isnan(gradpx).eq.1) lcomax=7.0d0
#elif machine == 7
if(ieee_is_nan(gradpx).eqv..true.) lcomax=7.0d0
#elif machine == 15
if(ieee_is_nan(gradpx).eqv..true.) lcomax=7.0d0
#elif machine == 16
if(ieee_is_nan(gradpx).eqv..true.) lcomax=7.0d0
#elif machine == 17
if(ieee_is_nan(gradpx).eqv..true.) lcomax=7.0d0
#else
if(isnan(gradpx).eqv..true.) lcomax=7.0d0
#endif

#if phi_flag == 1
#if machine == 4
if(isnan(phic(1,1,1,1)).eq.1) lcomax=7.0d0
#elif machine == 7
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#elif machine == 15
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#elif machine == 16
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#elif machine == 17
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#else
if(isnan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#endif

#if psi_flag == 1
#if machine == 4
if(isnan(psic_fg(1,1,1,1)).eq.1) lcomax=7.0d0
#elif machine == 7
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#elif machine == 15
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#elif machine == 16
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#elif machine == 17
if(ieee_is_nan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#else
if(isnan(phic(1,1,1,1)).eqv..true.) lcomax=7.0d0
#endif
#endif
#endif

call mpi_allreduce(lcomax,gcomax,1,mpi_double,mpi_max,flow_comm,ierr)

if(rank.eq.0) write(*,'(1x,a,es8.2)') 'check on Courant number : ',gcomax

#if machine == 4
if((gcomax.gt.co).or.(isnan(gcomax).eq.1))then
  if(rank.eq.0) write(*,'(1x,a,es8.2,a,f8.3)') 'Courant number exceeds maximum value : ',gcomax,'>',co
  call exit(0)
endif
#elif machine == 7
if((gcomax.gt.co).or.(ieee_is_nan(gcomax).eq.1))then
  if(rank.eq.0) write(*,'(1x,a,es8.2,a,f8.3)') 'Courant number exceeds maximum value : ',gcomax,'>',co
  call exit(0)
endif
#elif machine == 15
if((gcomax.gt.co).or.(ieee_is_nan(gcomax).eq.1))then
  if(rank.eq.0) write(*,'(1x,a,es8.2,a,f8.3)') 'Courant number exceeds maximum value : ',gcomax,'>',co
  call exit(0)
endif
#elif machine == 16
if((gcomax.gt.co).or.(ieee_is_nan(gcomax).eq.1))then
  if(rank.eq.0) write(*,'(1x,a,es8.2,a,f8.3)') 'Courant number exceeds maximum value : ',gcomax,'>',co
  call exit(0)
endif
#elif machine == 17
if((gcomax.gt.co).or.(ieee_is_nan(gcomax).eq.1))then
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
