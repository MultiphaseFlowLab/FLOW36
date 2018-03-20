subroutine xz2xy(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata

!type(mpi_request) :: req,reqs(2)
integer :: req,iadd


integer :: dims(2)
integer :: source,dest,disp,direction,numel,indy,indz,ry,rz
integer :: sendy,recvz
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(nsxx,nz,npyy,2),wa(nsxx,npzz,ny,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)

! Vesta needs mpi_async_protects_nonblocking to be logical
!onlyforvesta


ry=mod(ny,nzcpu)
rz=mod(nz,nzcpu)

allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))

direction=0

do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)

 numel=ngxx*ngyy*ngzz*2
 bufs=0.0d0*bufs

 if((floor(real(dest)/real(nycpu)).lt.ry).or.ry.eq.0)then
  indy=floor(real(dest)/real(nycpu))*ngyy
  sendy=ngyy
 else
  indy=(floor(real(dest)/real(nycpu))-ry)*(ngyy-1)+ry*ngyy
  sendy=ngyy-1
 endif

 bufs(1:nsxx,1:npzz,1:sendy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+sendy,1:2)

! isend + recv  (blocking recv needs no wait)
 call mpi_isend(bufs,numel,mpi_double_precision,dest,16,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,16,cart_comm,mpi_status_ignore,ierr)


 if((floor(real(source)/real(nycpu)).lt.rz).or.rz.eq.0)then
  indz=floor(real(source)/real(nycpu))*ngzz
  recvz=ngzz
 else
  indz=(floor(real(source)/real(nycpu))-rz)*(ngzz-1)+rz*ngzz
  recvz=ngzz-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)


! write(*,*) 'rank',rank,'to',dest,indy+1,indy+sendy,'from',source,indz+1,indz+recvz
 uc(1:nsxx,indz+1:indz+recvz,1:npyy,1:2)=bufr(1:nsxx,1:recvz,1:npyy,1:2)
enddo

! copy data in place (avoid communication rank n to rank n)
if((floor(real(rank)/real(nycpu)).lt.ry).or.ry.eq.0)then
 indy=floor(real(rank)/real(nycpu))*ngyy
else
 indy=(floor(real(rank)/real(nycpu))-ry)*(ngyy-1)+ry*ngyy
endif
if((floor(real(rank)/real(nycpu)).lt.rz).or.rz.eq.0)then
 indz=floor(real(rank)/real(nycpu))*ngzz
else
 indz=(floor(real(rank)/real(nycpu))-rz)*(ngzz-1)+rz*ngzz
endif
!write(*,*) 'rank',rank,'to',rank,indy+1,indy+npyy,'from',rank,indz+1,indz+npzz
uc(1:nsxx,indz+1:indz+npzz,1:npyy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)

deallocate(bufs)
deallocate(bufr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine xz2xy_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata
use dual_grid

!type(mpi_request) :: req,reqs(2)
integer :: req,iadd


integer :: dims(2)
integer :: source,dest,disp,direction,numel,indy,indz,ry,rz
integer :: sendy,recvz
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(nsxx,npsiz,npyy,2),wa(nsxx,npzz,npsiy,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)

! Vesta needs mpi_async_protects_nonblocking to be logical
!onlyforvesta


ry=mod(npsiy,nzcpu)
rz=mod(npsiz,nzcpu)

allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))

direction=0

do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)

 numel=ngxx*ngyy*ngzz*2
 bufs=0.0d0*bufs

 if((floor(real(dest)/real(nycpu)).lt.ry).or.ry.eq.0)then
  indy=floor(real(dest)/real(nycpu))*ngyy
  sendy=ngyy
 else
  indy=(floor(real(dest)/real(nycpu))-ry)*(ngyy-1)+ry*ngyy
  sendy=ngyy-1
 endif

 bufs(1:nsxx,1:npzz,1:sendy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+sendy,1:2)

! isend + recv  (blocking recv needs no wait)
 call mpi_isend(bufs,numel,mpi_double_precision,dest,16,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,16,cart_comm,mpi_status_ignore,ierr)


 if((floor(real(source)/real(nycpu)).lt.rz).or.rz.eq.0)then
  indz=floor(real(source)/real(nycpu))*ngzz
  recvz=ngzz
 else
  indz=(floor(real(source)/real(nycpu))-rz)*(ngzz-1)+rz*ngzz
  recvz=ngzz-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)


! write(*,*) 'rank',rank,'to',dest,indy+1,indy+sendy,'from',source,indz+1,indz+recvz
 uc(1:nsxx,indz+1:indz+recvz,1:npyy,1:2)=bufr(1:nsxx,1:recvz,1:npyy,1:2)
enddo

! copy data in place (avoid communication rank n to rank n)
if((floor(real(rank)/real(nycpu)).lt.ry).or.ry.eq.0)then
 indy=floor(real(rank)/real(nycpu))*ngyy
else
 indy=(floor(real(rank)/real(nycpu))-ry)*(ngyy-1)+ry*ngyy
endif
if((floor(real(rank)/real(nycpu)).lt.rz).or.rz.eq.0)then
 indz=floor(real(rank)/real(nycpu))*ngzz
else
 indz=(floor(real(rank)/real(nycpu))-rz)*(ngzz-1)+rz*ngzz
endif
!write(*,*) 'rank',rank,'to',rank,indy+1,indy+npyy,'from',rank,indz+1,indz+npzz
uc(1:nsxx,indz+1:indz+npzz,1:npyy,1:2)=wa(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)

deallocate(bufs)
deallocate(bufr)

return
end
