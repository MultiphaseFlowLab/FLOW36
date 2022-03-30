subroutine spectral_to_phys(uc,uout,aliasing)

use commondata
use par_size
use mpi
use fftx_bwd_module

integer :: dims(2)
integer :: rx,ry,rz
integer :: aliasing
double precision :: uc(spx,nz,spy,2),uout(nx,fpz,fpy)
double precision, allocatable :: u(:,:,:,:),wa(:,:,:,:)


dims(1)=nzcpu
dims(2)=nycpu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    idct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(u(spx,nz,spy,2))
u=uc
call dctz_bwd(u,u,spx,nz,spy,aliasing)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2)    change parallelization x-y to x-z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rx=mod(nx/2+1,nycpu)
ry=mod(ny,nzcpu)
rz=mod(nz,nzcpu)
npz=int((nz-rz)/nzcpu)

if(rz.eq.0)then
 ngz=npz
else
 ngz=npz+1
endif

if(floor(real(rank)/real(nycpu)).lt.rz)then
 npz=int((nz-rz)/nzcpu)+1
endif
npx=spx
npy=spy

if(rx.eq.0)then
 ngx=spx
else
 ngx=spx+1
endif

if(mod(rank,nycpu).lt.rx)then
 ngx=spx
endif

if(ry.eq.0)then
 ngy=spy
else
 ngy=spy+1
endif

if(floor(real(rank)/real(nycpu)).lt.ry)then
 ngy=spy
endif


#define nzcpu nnzcpu
#if nzcpu>1
!if(nzcpu.gt.1)then ! substituted with conditional compilation

 allocate(wa(spx,nz,spy,2))
 wa=u
 deallocate(u)
 allocate(u(spx,npz,ny,2))

 call xy2xz(wa,u,dims,ngx,npx,ngy,npy,ngz,npz)

 deallocate(wa)

!endif
#else
!if(rank.eq.0) write(*,*) 'skipping xy - xz'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    ifft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_bwd(u,u,spx,npz,ny,aliasing)
!call ffty_bwd(u,u,spx,npz,ny,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4)    change parallelization x-z to y-z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ry=mod(ny,nycpu)
npy=int((ny-ry)/nycpu)

if(ry.eq.0)then
 ngy=npy
else
 ngy=npy+1
endif

if(mod(rank,nycpu).lt.ry) then
 npy=int((ny-ry)/nycpu)+1
endif


#define nycpu nnycpu
#if nycpu>1
!if(nycpu.gt.1)then ! substituted with conditional compilation

 allocate(wa(spx,npz,ny,2))
 wa=u
 deallocate(u)
 allocate(u(nx/2+1,npz,npy,2))

 call xz2yz(wa,u,dims,ngx,npx,ngy,npy,ngz,npz)

 deallocate(wa)

!endif
#else
!if(rank.eq.0) write(*,*) 'skipping xz - yz'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    ifft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call fftx_bwd(u,uout,aliasing)


deallocate(u)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spectral_to_phys_fg(uc,uout,aliasing)

use commondata
use par_size
use dual_grid
use mpi
use fftx_bwd_module

integer :: dims(2)
integer :: rx,ry,rz
integer :: aliasing

double precision :: uc(spxpsi,npsiz,spypsi,2),uout(npsix,fpzpsi,fpypsi)
double precision, allocatable :: u(:,:,:,:),wa(:,:,:,:)


dims(1)=nzcpu
dims(2)=nycpu


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    idct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(u(spxpsi,npsiz,spypsi,2))
u=uc
call dctz_bwd_fg(u,u,spxpsi,npsiz,spypsi,aliasing)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2)    change parallelization x-y to x-z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rx=mod(npsix/2+1,nycpu)
ry=mod(npsiy,nzcpu)
rz=mod(npsiz,nzcpu)
npz=int((npsiz-rz)/nzcpu)

if(rz.eq.0)then
 ngz=npz
else
 ngz=npz+1
endif

if(floor(real(rank)/real(nycpu)).lt.rz)then
 npz=int((npsiz-rz)/nzcpu)+1
endif
npx=spxpsi
npy=spypsi

if(rx.eq.0)then
 ngx=spxpsi
else
 ngx=spxpsi+1
endif

if(mod(rank,nycpu).lt.rx)then
 ngx=spxpsi
endif

if(ry.eq.0)then
 ngy=spypsi
else
 ngy=spypsi+1
endif

if(floor(real(rank)/real(nycpu)).lt.ry)then
 ngy=spypsi
endif


#define nzcpu nnzcpu
#if nzcpu>1

 allocate(wa(spxpsi,npsiz,spypsi,2))
 wa=u
 deallocate(u)
 allocate(u(spxpsi,npz,npsiy,2))

 call xy2xz_fg(wa,u,dims,ngx,npx,ngy,npy,ngz,npz)

 deallocate(wa)

#else
!if(rank.eq.0) write(*,*) 'skipping xy - xz'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    ifft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_bwd_fg(u,u,spxpsi,npz,npsiy,aliasing)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4)    change parallelization x-z to y-z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ry=mod(npsiy,nycpu)
npy=int((npsiy-ry)/nycpu)

if(ry.eq.0)then
 ngy=npy
else
 ngy=npy+1
endif

if(mod(rank,nycpu).lt.ry) then
 npy=int((npsiy-ry)/nycpu)+1
endif


#define nycpu nnycpu
#if nycpu>1

 allocate(wa(spxpsi,npz,npsiy,2))
 wa=u
 deallocate(u)
 allocate(u(npsix/2+1,npz,npy,2))

 call xz2yz_fg(wa,u,dims,ngx,npx,ngy,npy,ngz,npz)

 deallocate(wa)

#else
!if(rank.eq.0) write(*,*) 'skipping xz - yz'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    ifft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call fftx_bwd_fg(u,uout,aliasing)

deallocate(u)


return
end
