subroutine wave_numbers

use commondata
use wavenumber
use sim_par
use phase_field
use dual_grid

integer :: i,j

#define match_dens matched_density
#define match_visc matched_viscosity

kx(1)=0.0d0
do i=2,nx/2+1
  kx(i)=dble(i-1)*2.0d0*pi/xl
enddo

ky(1)=0.0d0
do i=2,ny/2+1
  ky(ny-i+2)=-dble(i-1)*2.0d0*pi/yl
  ky(i)=dble(i-1)*2.0d0*pi/yl
enddo

do j=1,ny
  do i=1,nx/2+1
    k2(i,j)=kx(i)*kx(i)+ky(j)*ky(j)
  enddo
enddo

allocate(kxpsi(npsix/2+1))
allocate(kypsi(npsiy))
allocate(k2psi(npsix/2+1,npsiy))

kxpsi(1)=0.0d0
do i=2,npsix/2+1
  kxpsi(i)=dble(i-1)*2.0d0*pi/xl
enddo

kypsi(1)=0.0d0
do i=2,npsiy/2+1
  kypsi(npsiy-i+2)=-dble(i-1)*2.0d0*pi/yl
  kypsi(i)=dble(i-1)*2.0d0*pi/yl
enddo

do j=1,npsiy
  do i=1,npsix/2+1
    k2psi(i,j)=kxpsi(i)*kxpsi(i)+kypsi(j)*kypsi(j)
  enddo
enddo


#if match_dens == 2
gamma=dt/(2.0d0*re*rhor)
#else
gamma=dt/(2.0d0*re)
#endif

#if match_visc == 2
gamma=gamma*visr
#endif

return
end
