subroutine spectral_to_phys(uct,uout,aliasing)

use commondata

integer :: aliasing

double precision :: uct(nx/2+1,nz,ny,2),uout(nx,nz,ny)
double precision, allocatable :: ut(:,:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    idct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(ut(nx/2+1,nz,ny,2))
ut=uct
call dctz_bwd(ut,ut,nx/2+1,nz,ny,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    ifft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_bwd(ut,ut,nx/2+1,nz,ny,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    ifft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call fftx_bwd(ut,uout,nx,nz,ny,aliasing)

deallocate(ut)


return
end
