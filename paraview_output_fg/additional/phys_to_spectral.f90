subroutine phys_to_spectral(ut,uout,aliasing)

use commondata

integer :: aliasing

double precision :: ut(nx,nz,ny),uout(nx/2+1,nz,ny,2)
double precision, allocatable :: uct(:,:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(uct(nx/2+1,nz,ny,2))
call fftx_fwd(ut,uct,nx,nz,ny,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_fwd(uct,uct,nx/2+1,nz,ny,aliasing) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    dct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call dctz_fwd(uct,uout,nx/2+1,nz,ny,aliasing)
!uout=uc

deallocate(uct)

return
end
