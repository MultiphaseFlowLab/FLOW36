subroutine phys_to_spectral(u,uout,aliasing)

use commondata
use mpi

integer :: npx
integer :: aliasing

double precision :: u(nx,nz,ny),uout(nx/2+1,nz,ny,2)
double precision, allocatable :: uc(:,:,:,:)

npx=nx/2+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(uc(npx,nz,ny,2))
call fftx_fwd(u,uc,nx,nz,ny,aliasing)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_fwd(uc,uc,npx,nz,ny,aliasing)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    dct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call dctz_fwd(uc,uout,npx,nz,ny,aliasing)


deallocate(uc)

return
end
