subroutine xy2xz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use mpi
use commondata

integer :: req,iadd

integer :: dims(2)
integer :: source,dest,disp,direction,numel,indz,indy,rz,ry
integer :: sendz,recvy
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz

double precision :: uc(nsxx,npzz,ny,2),wa(nsxx,nz,npyy,2)
double precision, allocatable, asynchronous :: bufs(:,:,:,:),bufr(:,:,:,:)


direction=0

ry=mod(ny,nzcpu)
rz=mod(nz,nzcpu)

allocate(bufs(ngxx,ngzz,ngyy,2))
allocate(bufr(ngxx,ngzz,ngyy,2))


do disp=1,dims(direction+1)-1
 call mpi_cart_shift(cart_comm,direction,disp,source,dest,ierr)

 numel=ngxx*ngyy*ngzz*2
 bufs=0.0d0*bufs

 if((floor(real(dest)/real(nycpu)).lt.rz).or.rz.eq.0)then
  indz=floor(real(dest)/real(nycpu))*ngzz
  sendz=ngzz
 else
  indz=(floor(real(dest)/real(nycpu))-rz)*(ngzz-1)+rz*ngzz
  sendz=ngzz-1
 endif

 bufs(1:nsxx,1:sendz,1:npyy,1:2)=wa(1:nsxx,indz+1:indz+sendz,1:npyy,1:2)

! isend + recv
 call mpi_isend(bufs,numel,mpi_double_precision,dest,18,cart_comm,req,ierr)
 call mpi_recv(bufr,numel,mpi_double_precision,source,18,cart_comm,mpi_status_ignore,ierr)


 if((floor(real(source)/real(nycpu)).lt.ry).or.ry.eq.0)then
  indy=floor(real(source)/real(nycpu))*ngyy
  recvy=ngyy
 else
  indy=(floor(real(source)/real(nycpu))-ry)*(ngyy-1)+ry*ngyy
  recvy=ngyy-1
 endif

! for use mpi_f08
! if(.not.mpi_async_protects_nonblocking) call mpi_f_sync_reg(bufs,ierr)
 call mpi_wait(req,mpi_status_ignore,ierr)
 if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufs,iadd,ierr)

! write(*,*) 'rank',rank,'to',dest,indz+1,indz+sendz,'from',source,indy+1,indy+recvy
 uc(1:nsxx,1:npzz,indy+1:indy+recvy,1:2)=bufr(1:nsxx,1:npzz,1:recvy,1:2)
enddo


! copy data in place (avoid communication rank n to rank n)
if((floor(real(rank)/real(nycpu)).lt.rz).or.rz.eq.0)then
 indz=floor(real(rank)/real(nycpu))*ngzz
else
 indz=(floor(real(rank)/real(nycpu))-rz)*(ngzz-1)+rz*ngzz
endif
if((floor(real(rank)/real(nycpu)).lt.ry).or.ry.eq.0)then
 indy=floor(real(rank)/real(nycpu))*ngyy
else
 indy=(floor(real(rank)/real(nycpu))-ry)*(ngyy-1)+ry*ngyy
endif

!write(*,*) 'rank',rank,'to',rank,indz+1,indz+npzz,'from',rank,indy+1,indy+npyy
uc(1:nsxx,1:npzz,indy+1:indy+npyy,1:2)=wa(1:nsxx,indz+1:indz+npzz,1:npyy,1:2)








deallocate(bufs)
deallocate(bufr)


return
end
