subroutine phys_to_spectral(u,uout,aliasing)

use commondata
use par_size
use mpi


integer :: dims(2) !,coord(2)
integer :: rx,ry,rz,nsx,ngsx
integer :: ngx,ngy,ngz,npx,npy,npz
integer :: aliasing

!double precision :: stime,etime,dtime,mtime
double precision :: u(nx,fpz,fpy),uout(spx,nz,spy,2)
double precision, allocatable :: uc(:,:,:,:),wa(:,:,:,:)

!! just to check code
!if(rank.eq.0) then
! open(unit=22,file='./aux_matlab/before.dat',form='formatted',status='replace')
! write(22,'(f16.8)') u(:,1,1)
! close(22,status='keep')
!endif
!!!!!!!!!!!!!! end

ry=mod(ny,nycpu)
rz=mod(nz,nzcpu)

npx=nx/2+1

npy=int((ny-ry)/nycpu)
if(ry.eq.0)then
 ngy=npy
else
 ngy=npy+1
endif
if(mod(rank,nycpu).lt.ry)then
 npy=int((ny-ry)/nycpu)+1
endif

npz=int((nz-rz)/nzcpu)
if(rz.eq.0)then
 ngz=npz
else
 ngz=npz+1
endif

if(floor(real(rank)/real(nycpu)).lt.rz)then
 npz=int((nz-rz)/nzcpu)+1
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1)    fft x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(uc(npx,npz,npy,2))
call fftx_fwd(u,uc,nx,npz,npy,aliasing)
!call fftx_fwd(u,uc,nx,npz,npy,0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2)    change parallelization y-z to x-z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


dims(1)=nzcpu
dims(2)=nycpu

!stime=mpi_wtime()


rx=mod(nx/2+1,nycpu)
nsx=int((npx-rx)/nycpu)
if(rx.eq.0)then
 ngsx=nsx
else
 ngsx=nsx+1
endif

if(mod(rank,nycpu).lt.rx)then
 nsx=int((nx/2+1-rx)/nycpu)+1
endif
ngx=ngsx
npx=nsx


#define nycpu nnycpu
#if nycpu>1
!if(nycpu.gt.1)then ! substituted with conditional compilation

 allocate(wa(nx/2+1,npz,npy,2))
 wa=uc
 deallocate(uc)
 allocate(uc(nsx,npz,ny,2))


 call yz2xz(wa,uc,dims,ngx,npx,ngy,npy,ngz,npz)

 deallocate(wa)

!endif
#else
!if(rank.eq.0) write(*,*) 'skipping yz - xz'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3)    fft y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ffty_fwd(uc,uc,nsx,npz,ny,aliasing) 
!call ffty_fwd(uc,uc,nsx,npz,ny,0) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4)    change parallelization x-z to x-y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ry=mod(ny,nzcpu)
npy=int((ny-ry)/nzcpu)
if(ry.eq.0)then
 ngy=npy
else
 ngy=npy+1
endif

if(floor(real(rank)/real(nycpu)).lt.ry)then 
 npy=int((ny-ry)/nzcpu)+1
endif
rz=mod(nz,nzcpu)

#define nzcpu nnzcpu
#if nzcpu>1
!if(nzcpu.gt.1)then ! substituted with conditional compilation

 allocate(wa(nsx,npz,ny,2))
 wa=uc
 deallocate(uc)
 allocate(uc(nsx,nz,npy,2))

 call xz2xy(wa,uc,dims,ngx,npx,ngy,npy,ngz,npz)

 deallocate(wa)

!endif
#else
!if(rank.eq.0) write(*,*) 'skipping xz - xy'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 5)    dct z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call dctz_fwd(uc,uout,nsx,nz,npy,aliasing)
!uout=uc


deallocate(uc)


!etime=mpi_wtime()
!dtime=etime-stime

!call mpi_reduce(dtime,mtime,1,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
!if(rank.eq.0)write(*,'(1x,A,F12.6,A)') 'mean communication time=',mtime/dble(ntask),' s'


return
end
