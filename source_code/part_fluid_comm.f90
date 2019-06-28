subroutine get_velocity

use mpi
use commondata
use velocity
use particle
use comm_pattern

integer :: i,number
double precision, allocatable :: bufs(:,:,:),bufr(:,:,:,:)


if(rank.eq.leader) allocate(bufr(chunk_size(1),chunk_size(2),chunk_size(3),chunk_size(4)))

number=chunk_size(1)*chunk_size(2)*chunk_size(3)

! get velocity, only ranks in comm_comm are involved
! send u
if(rank.le.flow_comm_lim)then
 allocate(bufs(chunk_size(1),chunk_size(2),chunk_size(3)))
 bufs=0.0d0
 if(rank.lt.flow_comm_lim)then
  bufs(:,1:saved_size(rank+1,2),1:saved_size(rank+1,3))=u
 endif
 call mpi_gather(bufs,number,mpi_double_precision,bufr,number,mpi_double_precision,leader,comm_comm,ierr)
 deallocate(bufs)
endif

! synchronization barrier
if(rank.ge.leader)then
 if(rank.eq.leader)then
  ! write data
  do i=1,flow_comm_lim
   uf(:,address_start(i,2)+1:address_start(i,2)+saved_size(i,2),address_start(i,3)+1:address_start(i,3)+saved_size(i,3))= &
 &       bufr(:,1:saved_size(i,2),1:saved_size(i,3),i)
  enddo
 endif
 call mpi_win_fence(0,window_u,ierr)
endif


! send v
if(rank.le.flow_comm_lim)then
 allocate(bufs(chunk_size(1),chunk_size(2),chunk_size(3)))
 bufs=0.0d0
 if(rank.lt.flow_comm_lim)then
  bufs(:,1:saved_size(rank+1,2),1:saved_size(rank+1,3))=v
 endif
 call mpi_gather(bufs,number,mpi_double_precision,bufr,number,mpi_double_precision,leader,comm_comm,ierr)
 deallocate(bufs)
endif

! synchronization barrier
if(rank.ge.leader)then
 if(rank.eq.leader)then
  ! write data
  do i=1,flow_comm_lim
   vf(:,address_start(i,2)+1:address_start(i,2)+saved_size(i,2),address_start(i,3)+1:address_start(i,3)+saved_size(i,3))= &
 &       bufr(:,1:saved_size(i,2),1:saved_size(i,3),i)
  enddo
 endif
 call mpi_win_fence(0,window_v,ierr)
endif


! send w
if(rank.le.flow_comm_lim)then
 allocate(bufs(chunk_size(1),chunk_size(2),chunk_size(3)))
 bufs=0.0d0
 if(rank.lt.flow_comm_lim)then
  bufs(:,1:saved_size(rank+1,2),1:saved_size(rank+1,3))=w
 endif
 call mpi_gather(bufs,number,mpi_double_precision,bufr,number,mpi_double_precision,leader,comm_comm,ierr)
 deallocate(bufs)
endif

! synchronization barrier
if(rank.ge.leader)then
 if(rank.eq.leader)then
  ! write data
  do i=1,flow_comm_lim
   wf(:,address_start(i,2)+1:address_start(i,2)+saved_size(i,2),address_start(i,3)+1:address_start(i,3)+saved_size(i,3))= &
 &       bufr(:,1:saved_size(i,2),1:saved_size(i,3),i)
  enddo
  deallocate(bufr)
 endif
 call mpi_win_fence(0,window_w,ierr)
endif

! if(rank.eq.leader+1) write(*,*) maxval(uf),minval(uf),shape(uf)

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_2WCforces

use mpi
use commondata
use velocity
use particle
use comm_pattern

integer :: i,number
double precision, allocatable :: bufs(:,:,:,:),bufr(:,:,:)


number=chunk_size(1)*chunk_size(2)*chunk_size(3)

! fill in Fb_x
if(rank.eq.leader)then
 allocate(bufs(chunk_size(1),chunk_size(2),chunk_size(3),chunk_size(4)))
 bufs=0.0d0
 do i=1,flow_comm_lim
  bufs(:,1:saved_size(i,2),1:saved_size(i,3),i)= &
 & fb_x(:,address_start(i,2)+1:address_start(i,2)+saved_size(i,2),address_start(i,3)+1:address_start(i,3)+saved_size(i,3))
 enddo
endif

if(rank.le.flow_comm_lim)then
 allocate(bufr(chunk_size(1),chunk_size(2),chunk_size(3)))
 bufr=0.0d0
 call mpi_scatter(bufs,number,mpi_double_precision,bufr,number,mpi_double_precision,leader,comm_comm,ierr)
 ! write bufr to force matrix
 if(rank.lt.flow_comm_lim)then
  forx=bufr(:,1:saved_size(rank+1,2),1:saved_size(rank+1,3))
 endif
 deallocate(bufr)
endif

! fill in Fb_y
if(rank.eq.leader)then
 allocate(bufs(chunk_size(1),chunk_size(2),chunk_size(3),chunk_size(4)))
 bufs=0.0d0
 do i=1,flow_comm_lim
  bufs(:,1:saved_size(i,2),1:saved_size(i,3),i)= &
 & fb_y(:,address_start(i,2)+1:address_start(i,2)+saved_size(i,2),address_start(i,3)+1:address_start(i,3)+saved_size(i,3))
 enddo
endif

if(rank.le.flow_comm_lim)then
 allocate(bufr(chunk_size(1),chunk_size(2),chunk_size(3)))
 bufr=0.0d0
 call mpi_scatter(bufs,number,mpi_double_precision,bufr,number,mpi_double_precision,leader,comm_comm,ierr)
 ! write bufr to force matrix
 if(rank.lt.flow_comm_lim)then
  fory=bufr(:,1:saved_size(rank+1,2),1:saved_size(rank+1,3))
 endif
 deallocate(bufr)
endif

! fill in Fb_z
if(rank.eq.leader)then
 allocate(bufs(chunk_size(1),chunk_size(2),chunk_size(3),chunk_size(4)))
 bufs=0.0d0
 do i=1,flow_comm_lim
  bufs(:,1:saved_size(i,2),1:saved_size(i,3),i)= &
 & fb_z(:,address_start(i,2)+1:address_start(i,2)+saved_size(i,2),address_start(i,3)+1:address_start(i,3)+saved_size(i,3))
 enddo
endif

if(rank.le.flow_comm_lim)then
 allocate(bufr(chunk_size(1),chunk_size(2),chunk_size(3)))
 bufr=0.0d0
 call mpi_scatter(bufs,number,mpi_double_precision,bufr,number,mpi_double_precision,leader,comm_comm,ierr)
 ! write bufr to force matrix
 if(rank.lt.flow_comm_lim)then
  forz=bufr(:,1:saved_size(rank+1,2),1:saved_size(rank+1,3))
 endif
 deallocate(bufr)
endif


if(rank.eq.leader) deallocate(bufs)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_communication_pattern

use mpi
use commondata
use particle
use comm_pattern

integer :: rem
integer :: i,j,k,ii
integer :: savedim(nycpu,nzcpu,3)


! define sizes of chunks for communications
! communicated size
chunk_size(1)=nx

rem=mod(nz,nzcpu)
if(rem.gt.0)then
 chunk_size(2)=int((nz-rem)/nzcpu)+1
else
 chunk_size(2)=int(nz/nzcpu)
endif

rem=mod(ny,nycpu)
if(rem.gt.0)then
 chunk_size(3)=int((ny-rem)/nycpu)+1
else
 chunk_size(3)=int(ny/nycpu)
endif

if(flow_comm_lim.eq.leader)then
 ! nodes>1, allocate additional space for leader
 chunk_size(4)=flow_comm_lim+1
else
 ! nodes=1, leader is rank #0
 chunk_size(4)=flow_comm_lim
endif

! saved size, size of domain differs among tasks, not all that is communicated is kept
! rank corresponds to row number-1
allocate(saved_size(flow_comm_lim,3))
saved_size(:,1)=nx

do i=1,flow_comm_lim
 ! y dimension
 rem=mod(ny,nycpu)
 saved_size(i,3)=int((ny-rem)/nycpu)
 if(mod(i-1,nycpu).lt.rem) saved_size(i,3)=saved_size(i,3)+1
 ! z dimension
 rem=mod(nz,nzcpu)
 saved_size(i,2)=int((nz-rem)/nzcpu)
 if(floor(dble(i-1)/dble(nycpu)).lt.rem) saved_size(i,2)=saved_size(i,2)+1
enddo

! write(*,*) 'Chunk ',chunk_size
! write(*,*) 'Saved ',rank,saved_size(rank+1,:)

! define mapping from cartesian geometry to linear array
! ranks are defined as iy-1+nycpu*(iz-1); iy=1:nycpu: column number, iz=1:nzcpu: row number

savedim=reshape(saved_size,[nycpu,nzcpu,3])

allocate(address_start(flow_comm_lim,3))
address_start=0

do k=2,nzcpu
 do j=1,nycpu
  i=(j-1)+nycpu*(k-1) ! rank equivalent
  ii=(j-1)+nycpu*(k-2)
  address_start(i+1,2)=address_start(ii+1,2)+savedim(j,k-1,2)
 enddo
enddo

do k=1,nzcpu
 do j=2,nycpu
  i=(j-1)+nycpu*(k-1) ! rank equivalent
  ii=(j-2)+nycpu*(k-1)
  address_start(i+1,3)=address_start(ii+1,3)+savedim(j-1,k,3)
 enddo
enddo

! write(*,*) 'Address ',rank,address_start(rank+1,:)


return
end
