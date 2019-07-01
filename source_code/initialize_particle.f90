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
! each rank works on part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2)

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


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_particle

use mpi
use commondata
use sim_par
use particle

integer :: i

#define twowayc twowaycflag

! get fluid velocity for initial tracking
call get_velocity

! initialize particle position
if(in_cond_part_pos.eq.0)then
  if(rank.eq.0) write(*,*) 'Initializing random particle position'
  if(rank.eq.leader)then
   call random_number(xp)
   ! position in plus units, x,y,z
   xp(:,1)=xp(:,1)*xl*re
   xp(:,2)=xp(:,2)*yl*re
   xp(:,3)=(2.0d0*xp(:,3))*re   ! particles in [0,2*Re]
   ! xp(:,3)=(2.0d0*xp(:,3)-1.0d0)*re    ! particles in [-Re,+Re]
  endif
  ! synchronize shared memory window
  if(rank.ge.leader) call mpi_win_fence(0,window_xp,ierr)
else
  if(rank.eq.0) write(*,*) 'Dafuq? Check in_cond input value'
  stop
endif

! initialize particle velocity
if(in_cond_part_vel.eq.0)then
  if(rank.eq.0) write(*,*) 'Initializing zero particle velocity'
  if(rank.ge.leader)then
    up(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:)=0.0d0
    call mpi_win_fence(0,window_up,ierr)
  endif
elseif(in_cond_part_vel.eq.1)then
  if(rank.eq.0) write(*,*) 'Initializing fluid velocity at particle position'
  if(rank.ge.leader)then
   do i=part_index(rank_loc+1,1)+1,part_index(rank_loc+1,2)
    call lagran4(xp(i,:),up(i,:))
   enddo
   call mpi_win_fence(0,window_up,ierr)
  endif
else
  if(rank.eq.0) write(*,*) 'Dafuq? Check in_cond input value'
  stop
endif


#if twowayc == 1
call get_2WCforces
#endif

! ! debug only
! if(rank_loc.eq.2)then
!  open(456,file='./results/part.dat',status='new',form='formatted')
!  do i=1,part_number
!   write(456,'(3(2x,f12.4),4x,3(2x,f12.4))') xp(i,:),up(i,:)
!  enddo
!  close(456,status='keep')
! endif

return
end
