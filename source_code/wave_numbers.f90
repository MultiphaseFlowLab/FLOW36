subroutine wave_numbers

use commondata
use wavenumber
use sim_par
use phase_field

integer :: i,j

#define match_dens matched_density

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

#if match_dens == 2
gamma=dt/(2.0d0*re*rhor)
#else
gamma=dt/(2.0d0*re)
#endif

return
end
