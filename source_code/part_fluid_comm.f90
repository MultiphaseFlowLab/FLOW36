subroutine get_velocity

use mpi
use commondata
use particle


! get velocity, only ranks in comm_comm are involved
if(leader.eq.0)then
 ! everyone execute: comm_comm=mpi_comm_world


else
 ! only ranks in comm_comm execute
 if(rank.le.leader)then

 endif
endif


! synchronization barrier
if(rank.ge.leader)then
 call mpi_win_fence(0,window_u,ierr)
 call mpi_win_fence(0,window_v,ierr)
 call mpi_win_fence(0,window_w,ierr)
endif

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_2WCforces




return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_communication_pattern

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
