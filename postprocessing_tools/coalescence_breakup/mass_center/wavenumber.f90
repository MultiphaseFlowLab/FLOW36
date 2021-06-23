subroutine wavenumber

use commondata
use wavenumbers

allocate(kx(nx/2+1))
allocate(ky(ny))
allocate(k2(nx/2+1,ny))


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


return
end
