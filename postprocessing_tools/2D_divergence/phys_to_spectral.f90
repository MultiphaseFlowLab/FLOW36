subroutine phys_to_spectral(var,uout,aliasing)

use commondata

integer :: aliasing

double precision :: var(nx,nz,ny)
double precision :: uout(nx/2+1,nz,ny,2)
double precision, allocatable :: varc(:,:,:,:)

integer :: nfx,nfy,nfz

nfx=nx
nfy=ny
nfz=nz

npx=nfx/2+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(varc(npx,nfz,nfy,2))
call fftx_fwd(var,varc,nfx,nfz,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_fwd(varc,varc,npx,nfz,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    dct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call dctz_fwd(varc,uout,npx,nfz,nfy,aliasing)


deallocate(varc)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phys_to_spectral_fg(var,uout,aliasing)

use commondata

integer :: aliasing

double precision :: var(expx*nx,expz*(nz-1)+1,expy*ny)
double precision :: uout(expx*nx/2+1,expz*(nz-1)+1,expy*ny,2)
double precision, allocatable :: varc(:,:,:,:)

integer :: nfx,nfy,nfz

nfx=expx*nx
nfy=expy*ny
nfz=expz*(nz-1)+1

npx=nfx/2+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(varc(npx,nfz,nfy,2))
call fftx_fwd_fg(var,varc,nfx,nfz,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_fwd_fg(varc,varc,npx,nfz,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    dct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call dctz_fwd_fg(varc,uout,npx,nfz,nfy,aliasing)


deallocate(varc)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phys_to_spectral_2D(var,uout,aliasing)

use commondata

integer :: aliasing

double precision :: var(nx,ny)
double precision :: uout(nx/2+1,ny,2)
double precision, allocatable :: varc(:,:,:)

integer :: nfx,nfy

nfx=nx
nfy=ny

npx=nfx/2+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(varc(npx,nfy,2))
call fftx_fwd_2D(var,varc,nfx,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call ffty_fwd_2D(varc,uout,npx,nfy,aliasing)

deallocate(varc)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
