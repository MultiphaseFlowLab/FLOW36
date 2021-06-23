subroutine generate_grid

use commondata
use grid
use sim_par

double precision :: dx,dy

integer :: i,j,k

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)

x(1)=0.0d0
y(1)=0.0d0
do i=2,nx
  x(i)=x(i-1)+dx
enddo
do j=2,ny
  y(j)=y(j-1)+dy
enddo

do k =1,nz
  z(k)=dcos(((k-1)*pi)/(nz-1))
enddo

return
end
