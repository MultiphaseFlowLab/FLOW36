subroutine yz2xz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata

! for use mpi_f08
!type(mpi_request) :: req,reqs(2)
integer :: req,iadd  !,reqs(2)

integer :: dims(2)
integer :: source,dest,disp,direction,numel,indx,indy,rx,ry
integer :: sendx,recvy
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(nsxx,npzz,ny,2),wa(nx/2+1,npzz,npyy,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)

! Vesta needs mpi_async_protects_nonblocking to be logical
!onlyforvesta


! 0: z direction, 1: y direction (fortran columnwise order)
direction=1

rx=mod(nx/2+1,nycpu)
ry=mod(ny,nycpu)


allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))


do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)


! numel may be optimized, but now bufs and bufr are allocated just once
 numel=ngxx*ngyy*ngzz*2
 bufs=0.0d0*bufs

 if((mod(dest,nycpu).lt.rx).or.rx.eq.0)then
  indx=mod(dest,nycpu)*ngxx
  sendx=ngxx
 else
  indx=(mod(dest,nycpu)-rx)*(ngxx-1)+rx*ngxx
  sendx=ngxx-1
 endif


! do j=1,npy
!  do k=1,npz
!   do i=1,npx
!    bufs(i,k,j)=u(indx+i,k,j)
!   enddo
!  enddo
! enddo
 bufs(1:sendx,1:npzz,1:npyy,1:2)=wa(indx+1:indx+sendx,1:npzz,1:npyy,1:2)


! MPI communication among neighbours found, possible communications routines

! isend + recv
 call mpi_isend(bufs,numel,mpi_double_precision,dest,15,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,15,cart_comm,mpi_status_ignore,ierr)

!! isend + irecv + waitall
! call mpi_isend(bufs,numel,mpi_double_precision,dest,15,cart_comm,reqs(1),ierr)
! call mpi_irecv(bufr,numel,mpi_double_precision,source,15,cart_comm,reqs(2),ierr)
! call mpi_waitall(2,reqs,mpi_statuses_ignore,ierr)
! IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
!!call mpi_wait(req,mpi_status_ignore,ierr) ! mettere riga ... e buffer asincrono?

!! irecv + send
! call mpi_irecv(bufr,numel,mpi_double_precision,source,15,cart_comm,req,ierr)
! call mpi_send(bufs,numel,mpi_double_precision,dest,15,cart_comm,ierr)
! IF (.NOT.MPI_ASYNC_PROTECTS_NONBLOCKING) CALL MPI_F_sync_reg(snd_buf)
! call mpi_wait(req,mpi_status_ignore,ierr) ! mettere riga ... e buffer asincrono?



 if((mod(source,nycpu).lt.ry).or.ry.eq.0)then
  indy=mod(source,nycpu)*ngyy
  recvy=ngyy
 else
  indy=(mod(source,nycpu)-ry)*(ngyy-1)+ry*ngyy
  recvy=ngyy-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)

! write(*,*) 'rank',rank,'to',dest,indx+1,indx+sendx,'from',source,indy+1,indy+recvy
 uc(1:nsxx,1:npzz,indy+1:indy+recvy,1:2)=bufr(1:nsxx,1:npzz,1:recvy,1:2)
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
!write(*,*) 'rank',rank,'to',rank,indx+1,indx+nsxx,'from',rank,indy+1,indy+npyy
uc(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)=wa(indx+1:indx+nsxx,1:npzz,1:npyy,1:2)

deallocate(bufs)
deallocate(bufr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine yz2xz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata
use dual_grid

! for use mpi_f08
!type(mpi_request) :: req,reqs(2)
integer :: req,iadd  !,reqs(2)

integer :: dims(2)
integer :: source,dest,disp,direction,numel,indx,indy,rx,ry
integer :: sendx,recvy
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(nsxx,npzz,npsiy,2),wa(npsix/2+1,npzz,npyy,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)

! Vesta needs mpi_async_protects_nonblocking to be logical
!onlyforvesta


! 0: z direction, 1: y direction (fortran columnwise order)
direction=1

rx=mod(npsix/2+1,nycpu)
ry=mod(npsiy,nycpu)


allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))


do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)


! numel may be optimized, but now bufs and bufr are allocated just once
 numel=ngxx*ngyy*ngzz*2
 bufs=0.0d0*bufs

 if((mod(dest,nycpu).lt.rx).or.rx.eq.0)then
  indx=mod(dest,nycpu)*ngxx
  sendx=ngxx
 else
  indx=(mod(dest,nycpu)-rx)*(ngxx-1)+rx*ngxx
  sendx=ngxx-1
 endif

 bufs(1:sendx,1:npzz,1:npyy,1:2)=wa(indx+1:indx+sendx,1:npzz,1:npyy,1:2)


! MPI communication among neighbours found, possible communications routines

! isend + recv
 call mpi_isend(bufs,numel,mpi_double_precision,dest,15,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,15,cart_comm,mpi_status_ignore,ierr)

 if((mod(source,nycpu).lt.ry).or.ry.eq.0)then
  indy=mod(source,nycpu)*ngyy
  recvy=ngyy
 else
  indy=(mod(source,nycpu)-ry)*(ngyy-1)+ry*ngyy
  recvy=ngyy-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)

! write(*,*) 'rank',rank,'to',dest,indx+1,indx+sendx,'from',source,indy+1,indy+recvy
 uc(1:nsxx,1:npzz,indy+1:indy+recvy,1:2)=bufr(1:nsxx,1:npzz,1:recvy,1:2)
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
!write(*,*) 'rank',rank,'to',rank,indx+1,indx+nsxx,'from',rank,indy+1,indy+npyy
uc(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)=wa(indx+1:indx+nsxx,1:npzz,1:npyy,1:2)

deallocate(bufs)
deallocate(bufr)

return
end
