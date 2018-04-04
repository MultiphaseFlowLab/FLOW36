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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spectral_to_phys_fg(uct,uout,aliasing)

use commondata

integer :: aliasing
integer :: nfx,nfy,nfz

double precision :: uct((exp_x*nx)/2+1,exp_z*(nz-1)+1,exp_y*ny,2)
double precision :: uout(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny)
double precision, allocatable :: ut(:,:,:,:)

nfx=exp_x*nx
nfy=exp_y*ny
nfz=exp_z*(nz-1)+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    idct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(ut(nfx/2+1,nfz,nfy,2))
ut=uct
call dctz_bwd_fg(ut,ut,nfx/2+1,nfz,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    ifft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_bwd_fg(ut,ut,nfx/2+1,nfz,nfy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    ifft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call fftx_bwd_fg(ut,uout,nfx,nfz,nfy,aliasing)

deallocate(ut)


return
end
