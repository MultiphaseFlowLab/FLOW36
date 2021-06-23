subroutine wave_numbers

use commondata
use wavenumber
use sim_par
use phase_field

integer :: i,j

allocate(kx(nx/2+1))
allocate(ky(ny))

kx(1)=0.0d0
do i=2,nx/2+1
  kx(i)=dble(i-1)*2.0d0*pi/xl
enddo

ky(1)=0.0d0
do j=2,ny/2+1
  ky(ny-j+2)=-dble(j-1)*2.0d0*pi/yl
  ky(j)=dble(j-1)*2.0d0*pi/yl
enddo

return
end
