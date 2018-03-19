subroutine spectral_to_phys(uc,uout,aliasing)

use commondata
use par_size
use mpi

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


call fftx_bwd(u,uout,nx,fpz,fpy,aliasing)
!call fftx_bwd(u,uout,nx,fpz,fpy,0)

deallocate(u)


return
end
