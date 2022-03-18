program main

use commondata
use fields
use wavenumber

integer :: i,j

call read_input

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(u(nx,nz,ny))
allocate(v(nx,nz,ny))
allocate(w(nx,nz,ny))
allocate(phi(nx,nz,ny))
if(phi_flag.eq.0) phi=-1.0d0
allocate(kx(nx/2+1))
allocate(ky(ny))
allocate(k2(nx/2+1,ny))

! wavenumbers
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

! plan creation
call create_plan
call create_plan_2D

call read_grid

do i=nstart,nend,delta
 call divergence2D(i)
enddo
! run last step if not already run
if(mod(nend,delta).ne.0) call divergence2D(nend)



call destroy_plan

deallocate(x,y,z)
deallocate(u,v,w)
deallocate(phi)
deallocate(kx,ky,k2)


return
end program main
