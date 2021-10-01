subroutine phys_to_spectral(var,uout,aliasing)

use commondata


integer :: spx
integer :: aliasing

double precision :: uout(nx/2+1,nz,ny,2),var(nx,nz,ny)
double precision, allocatable :: varc(:,:,:,:)

spx=nx/2+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(varc(spx,nz,ny,2))
call fftx_fwd(var,varc,nx,nz,ny,aliasing)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call ffty_fwd(varc,varc,spx,nz,ny,aliasing)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    dct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call dctz_fwd(varc,uout,spx,nz,ny,aliasing)
! !uout=uc


deallocate(varc)


return
end
