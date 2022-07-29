subroutine xz2yz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata

! for use mpi_f08
!type(mpi_request) :: req,reqs(2)
integer :: req
integer :: iadd

integer :: dims(2)
integer :: source,dest,disp,direction,numel,indx,indy,rx,ry
integer :: sendy,recvx
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(nx/2+1,npzz,npyy,2),wa(nsxx,npzz,ny,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)

! Vesta needs mpi_async_protects_nonblocking to be logical
!onlyforvesta


direction=1


ry=mod(ny,nycpu)
rx=mod(nx/2+1,nycpu)

allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))

do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)

 numel=ngxx*ngyy*ngzz*2
 !$acc data copy(bufs)
 !$acc kernels
 bufs=0.0d0*bufs
 !$acc end kernels
 !$acc end data

 if((mod(dest,nycpu).lt.ry).or.ry.eq.0)then
  indy=mod(dest,nycpu)*ngyy
  sendy=ngyy
 else
  indy=(mod(dest,nycpu)-ry)*(ngyy-1)+ry*ngyy
  sendy=ngyy-1
 endif

 !$acc data copyin(wa) copyout(bufs)
 !$acc kernels
 bufs(1:nsxx,1:npzz,1:sendy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+sendy,1:2)
 !$acc end kernels
 !$acc end data

 ! isend + recv
 !CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
 !$acc host_data use_device(bufs,bufr)
 call mpi_isend(bufs,numel,mpi_double_precision,dest,17,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,17,cart_comm,mpi_status_ignore,ierr)
 !$acc end host_data

 if((mod(source,nycpu).lt.rx).or.rx.eq.0)then
  indx=mod(source,nycpu)*ngxx
  recvx=ngxx
 else
  indx=(mod(source,nycpu)-rx)*(ngxx-1)+rx*ngxx
  recvx=ngxx-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)

! write(*,*) 'rank',rank,'to',dest,indy+1,indy+sendy,'from',source,indx+1,indx+recvx
 !$acc data copyin(bufr) copyout(uc)
 !$acc kernels
 uc(indx+1:indx+recvx,1:npzz,1:npyy,1:2)=bufr(1:recvx,1:npzz,1:npyy,1:2)
 !$acc end kernels
 !$acc end data
enddo

! copy data in place (avoid communication rank n to rank n)
if((mod(rank,nycpu).lt.rx).or.rx.eq.0)then
 indx=mod(rank,nycpu)*ngxx
else
 indx=(mod(rank,nycpu)-rx)*(ngxx-1)+rx*ngxx
endif
if((mod(rank,nycpu).lt.ry).or.ry.eq.0)then
 indy=mod(rank,nycpu)*ngyy
else
 indy=(mod(rank,nycpu)-ry)*(ngyy-1)+ry*ngyy
endif
!write(*,*) 'rank',rank,'to',rank,indy+1,indy+npyy,'from',rank,indx+1,indx+nsxx
!$acc data copyin(wa) copyout(uc)
!$acc kernels
uc(indx+1:indx+nsxx,1:npzz,1:npyy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)
!$acc end kernels
!$acc end data

deallocate(bufs)
deallocate(bufr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine xz2yz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata
use dual_grid

! for use mpi_f08
!type(mpi_request) :: req,reqs(2)
integer :: req
integer :: iadd

integer :: dims(2)
integer :: source,dest,disp,direction,numel,indx,indy,rx,ry
integer :: sendy,recvx
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(npsix/2+1,npzz,npyy,2),wa(nsxx,npzz,npsiy,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)

! Vesta needs mpi_async_protects_nonblocking to be logical
!onlyforvesta

direction=1

ry=mod(npsiy,nycpu)
rx=mod(npsix/2+1,nycpu)

allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))

do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)

 numel=ngxx*ngyy*ngzz*2
 !$acc data copy(bufs)
 !$acc kernels
 bufs=0.0d0*bufs
 !$acc end kernels
 !$acc end data

 if((mod(dest,nycpu).lt.ry).or.ry.eq.0)then
  indy=mod(dest,nycpu)*ngyy
  sendy=ngyy
 else
  indy=(mod(dest,nycpu)-ry)*(ngyy-1)+ry*ngyy
  sendy=ngyy-1
 endif
 !$acc data copyin(wa) copyout(bufs)
 !$acc kernels
 bufs(1:nsxx,1:npzz,1:sendy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+sendy,1:2)
 !$acc end kernels
 !$acc end data

 ! isend + recv
 !CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
 !$acc host_data use_device(bufs,bufr)
 call mpi_isend(bufs,numel,mpi_double_precision,dest,17,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,17,cart_comm,mpi_status_ignore,ierr)
 !$acc end host_data

 if((mod(source,nycpu).lt.rx).or.rx.eq.0)then
  indx=mod(source,nycpu)*ngxx
  recvx=ngxx
 else
  indx=(mod(source,nycpu)-rx)*(ngxx-1)+rx*ngxx
  recvx=ngxx-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)

 ! write(*,*) 'rank',rank,'to',dest,indy+1,indy+sendy,'from',source,indx+1,indx+recvx
 !$acc data copyin(bufr) copyout(uc)
 !$acc kernels
 uc(indx+1:indx+recvx,1:npzz,1:npyy,1:2)=bufr(1:recvx,1:npzz,1:npyy,1:2)
 !$acc end kernels
 !$acc end data
enddo

! copy data in place (avoid communication rank n to rank n)
if((mod(rank,nycpu).lt.rx).or.rx.eq.0)then
 indx=mod(rank,nycpu)*ngxx
else
 indx=(mod(rank,nycpu)-rx)*(ngxx-1)+rx*ngxx
endif
if((mod(rank,nycpu).lt.ry).or.ry.eq.0)then
 indy=mod(rank,nycpu)*ngyy
else
 indy=(mod(rank,nycpu)-ry)*(ngyy-1)+ry*ngyy
endif
!write(*,*) 'rank',rank,'to',rank,indy+1,indy+npyy,'from',rank,indx+1,indx+nsxx
!$acc data copyin(wa) copyout(uc)
!$acc kernels
uc(indx+1:indx+nsxx,1:npzz,1:npyy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)
!$acc end kernels
!$acc end data

deallocate(bufs)
deallocate(bufr)

return
end
