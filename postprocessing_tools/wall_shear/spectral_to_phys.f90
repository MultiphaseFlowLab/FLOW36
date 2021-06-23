subroutine spectral_to_phys(uc,uout,aliasing)

use commondata
use mpi

integer :: aliasing,spx

double precision :: uc(nx/2+1,nz,ny,2),uout(nx,nz,ny)
double precision, allocatable :: u(:,:,:,:)

spx=nx/2+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    idct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(u(spx,nz,ny,2))
u=uc
call dctz_bwd(u,u,spx,nz,ny,aliasing)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    ifft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_bwd(u,u,spx,nz,ny,aliasing)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    ifft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call fftx_bwd(u,uout,nx,nz,ny,aliasing)


deallocate(u)


return
end
