subroutine allocate_particle

use, intrinsic :: iso_c_binding
use mpi
use commondata
use particle

type(c_ptr) :: baseptr
integer(kind=mpi_address_kind) :: varsize,lb
integer :: dispunit
integer :: number,delta,range
integer :: shape3(3),shape2(2)
integer :: i


if(rank.eq.leader)then
 ! only leader allocates memory, others get zero memory
 call mpi_type_get_extent(mpi_double_precision,lb,varsize,ierr)
 dispunit=int(varsize)
elseif(rank.gt.leader)then
 varsize=0
 dispunit=1
endif

! open share memory region for fluid flow, force feedback and particle data
shape2=[part_number,3]
shape3=[nx,nz,ny]

! allocate uf
number=nx*ny*nz
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_u,ierr)
! get location of memory segment
if(rank.gt.leader)then
  ! leader (global rank) is rank 0 in part_comm
  call mpi_win_shared_query(window_u,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,uf,shape3)

! allocate vf
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_v,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_v,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,vf,shape3)

! allocate wf
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_w,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_w,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,wf,shape3)

! allocate F_x
number=nx*ny*nz
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_fx,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_fx,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,fb_x,shape3)

! allocate F_y
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_fy,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_fy,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,fb_y,shape3)

! allocate F_z
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_fz,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_fz,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,fb_z,shape3)


! allocate particle data
number=part_number*3
! allocate xp
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_xp,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_xp,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,xp,shape2)

! allocate up
call mpi_win_allocate_shared(number*varsize,dispunit,mpi_info_null,part_comm,baseptr,window_up,ierr)
! get location of memory segment
if(rank.gt.leader)then
  call mpi_win_shared_query(window_up,0,number*varsize,dispunit,baseptr,ierr)
endif
! associate C pointer to Fortran pointer
call c_f_pointer(baseptr,up,shape2)


! define ranges of each part_comm processor within xp and up (start and end indexes)
allocate(part_index(ntask_sh,2))

! part index: from start_index+1 to start_index+delta
! -first column: start_index
! -second column: range of particles for each core (delta)




delta=mod(part_number,ntask_sh)
range=int((part_number-delta)/ntask_sh)
part_index(1,1)=0
part_index(:,2)=range
do i=1,ntask_sh
  if(i-1.lt.delta) part_index(i,2)=part_index(i,2)+1
enddo
do i=2,ntask_sh
  part_index(i,1)=part_index(i-1,1)+part_index(i-1,2)
enddo


! is it ok to use multiple shared windows? could it happen that they overwrite each other?
! can they correspond to the same physical address?
! if(rank.eq.leader) uf=-4.0d0
! if(rank.eq.leader+1) vf=-3.0d0
!
! ! synchronization
! call mpi_win_fence(0,window_u,ierr)
! call mpi_win_fence(0,window_v,ierr)
!
! if(rank.eq.leader) write(*,*) maxval(uf),minval(uf),maxval(vf),minval(vf)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_particle

use mpi
use commondata
use particle


! initialize particle position
if(in_cond_part_pos.eq.0)then
  if(rank.eq.0) write(*,*) 'Initialize random particle position'

else
  if(rank.eq.0) write(*,*) 'Dafuq? Check in_cond input value'
  stop
endif

! initialize particle velocity
if(in_cond_part_vel.eq.0)then
  if(rank.eq.0) write(*,*) 'Initialize zero particle velocity'

elseif(in_cond_part_vel.eq.1)then
  if(rank.eq.0) write(*,*) 'Initialize fluid velocity at particle position'

else
  if(rank.eq.0) write(*,*) 'Dafuq? Check in_cond input value'
  stop
endif


return
end
