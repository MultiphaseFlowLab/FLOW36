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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phys_to_spectral_fg(var,uout,aliasing)

use commondata

integer :: aliasing

double precision :: var(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny)
double precision :: uout(exp_x*nx/2+1,exp_z*(nz-1)+1,exp_y*ny,2)
double precision, allocatable :: varc(:,:,:,:)

integer :: nfx,nfy,nfz

nfx=exp_x*nx
nfy=exp_y*ny
nfz=exp_z*(nz-1)+1

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
