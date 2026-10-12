subroutine xy2xz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,ny,2),wa(nsxx,nz,npyy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nzcpu,8))
call xy2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,ny,nz)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine xy2xz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use dual_grid
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,npsiy,2),wa(nsxx,npsiz,npyy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nzcpu,8))
call xy2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,npsiy,npsiz)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! x-y to x-z pencils (direction 0, nzcpu ranks with the same y coordinate), gny x gnz grid
! wa: z complete, y split -> uc: y complete, z split
! Buffers: a2a_send/a2a_recv of module a2a_buffers, device memory allocated once (1D,
! block j at offset j*numel, ordered as (ngxx,ngzz,ngyy,2)); device (not managed) memory so
! that CUDA-aware MPI uses GPUDirect RDMA across nodes instead of staging through the host.
! All the blocks are packed with one kernel, exchanged with non-blocking messages on cart_comm_dir(0)
! and unpacked with one kernel. Blocks have the padded size ngxx*ngzz*ngyy*2 as in the previous
! pairwise exchange; block j goes to / comes from the rank with z coordinate j.
subroutine xy2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,gny,gnz)

use mpi
use commondata
use a2a_buffers

integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz,gny,gnz
integer :: me,nreq,req(2*nzcpu)
integer :: i,ky,kz,c,j,ry,rz,numel,i0,cnt
double precision :: uc(nsxx,npzz,gny,2),wa(nsxx,gnz,npyy,2)
integer(kind=8) :: nb,ib

ry=mod(gny,nzcpu)
rz=mod(gnz,nzcpu)

numel=ngxx*ngyy*ngzz*2
nb=int(numel,8)

! pack: block j holds the z slab of rank j (offset i0, cnt points), all local y
!$acc parallel loop collapse(5) private(i0,cnt,ib)
do j=0,nzcpu-1
 do c=1,2
  do ky=1,npyy
   do kz=1,ngzz
    do i=1,nsxx
     ib=int(j,8)*nb+int(i+ngxx*((kz-1)+ngzz*((ky-1)+ngyy*(c-1))),8)
     if(j.lt.rz .or. rz.eq.0)then
      i0=j*ngzz
      cnt=ngzz
     else
      i0=(j-rz)*(ngzz-1)+rz*ngzz
      cnt=ngzz-1
     endif
     if(kz.le.cnt) a2a_send(ib)=wa(i,i0+kz,ky,c)
    enddo
   enddo
  enddo
 enddo
enddo

! exchange with all the other ranks of cart_comm_dir(0) at once: all the receives and sends
! are posted together (no sequential rounds); the own block (j=me) is not sent, the unpack
! kernel takes it directly from a2a_send. Point-to-point instead of mpi_alltoall: the collective
! would copy the own block (and possibly stage blocks) on the host with managed memory.
call mpi_comm_rank(cart_comm_dir(0),me,ierr)
nreq=0
!CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
do j=0,nzcpu-1
 if(j.ne.me)then
  nreq=nreq+1
  call mpi_irecv(a2a_recv(1+int(j,8)*nb),numel,mpi_double_precision,j,0,cart_comm_dir(0),req(nreq),ierr)
 endif
enddo
do j=0,nzcpu-1
 if(j.ne.me)then
  nreq=nreq+1
  call mpi_isend(a2a_send(1+int(j,8)*nb),numel,mpi_double_precision,j,0,cart_comm_dir(0),req(nreq),ierr)
 endif
enddo
call mpi_waitall(nreq,req,mpi_statuses_ignore,ierr)

! unpack: block j holds the local z slab with the y points of rank j (offset i0, cnt points)
!$acc parallel loop collapse(5) private(i0,cnt,ib)
do j=0,nzcpu-1
 do c=1,2
  do ky=1,ngyy
   do kz=1,npzz
    do i=1,nsxx
     ib=int(j,8)*nb+int(i+ngxx*((kz-1)+ngzz*((ky-1)+ngyy*(c-1))),8)
     if(j.lt.ry .or. ry.eq.0)then
      i0=j*ngyy
      cnt=ngyy
     else
      i0=(j-ry)*(ngyy-1)+ry*ngyy
      cnt=ngyy-1
     endif
     if(ky.le.cnt)then
      if(j.eq.me)then
       uc(i,kz,i0+ky,c)=a2a_send(ib)
      else
       uc(i,kz,i0+ky,c)=a2a_recv(ib)
      endif
     endif
    enddo
   enddo
  enddo
 enddo
enddo

return
end
