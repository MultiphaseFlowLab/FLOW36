subroutine split_comm

use mpi
use commondata

integer :: shmem_comm
integer :: group_w,group_l
integer :: check_task
integer, allocatable, dimension(:) :: list_rank

integer :: i

#define particles particlecompflag

! query number of MPI processes, assign rank number (global rank and number)
call mpi_comm_size(mpi_comm_world,ntask_gl,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

#if particles == 1
  ! divide comm_world in shared memory commuicators
  call mpi_comm_split_type(mpi_comm_world,mpi_comm_type_shared,rank,mpi_info_null,shmem_comm,ierr)
  ! get size of shared memory communicator
  call mpi_comm_size(shmem_comm,ntask_sh,ierr)
  ! free communicator, will not be used
  call mpi_comm_free(shmem_comm,ierr)

  nodes=0
  do while(nodes*ntask_sh.lt.ntask_gl)
    nodes=nodes+1
  enddo

  ! create world group
  call mpi_comm_group(mpi_comm_world,group_w,ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! create particle communicator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! global rank of leader of part_comm
  leader=(nodes-1)*ntask_sh

  allocate(list_rank(ntask_sh))
  list_rank(1)=leader
  do i=2,ntask_sh
    list_rank(i)=list_rank(i-1)+1
  enddo
  call mpi_group_incl(group_w,ntask_sh,list_rank,group_l,ierr)
  deallocate(list_rank)

  call mpi_comm_create(mpi_comm_world,group_l,part_comm,ierr)

  call mpi_group_free(group_l,ierr)

  ! check that part_comm is a shared memory communicator
  if(rank.ge.leader)then
    call mpi_comm_split_type(part_comm,mpi_comm_type_shared,rank,mpi_info_null,shmem_comm,ierr)
    call mpi_comm_size(shmem_comm,check_task,ierr)
    call mpi_comm_rank(part_comm,rank_loc,ierr)
    if(check_task.ne.ntask_sh)then
      write(*,*) 'Particle communicator is not a shared memory communicator, stopping'
      call exit(0)
    endif
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! create flow communicator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(leader.gt.0)then
    ! at least 2 shared memory partitions (nodes>1)
    allocate(list_rank(leader))
    list_rank(1)=0
    do i=2,leader
      list_rank(i)=list_rank(i-1)+1
    enddo
    call mpi_group_incl(group_w,leader,list_rank,group_l,ierr)
    deallocate(list_rank)
    call mpi_comm_create(mpi_comm_world,group_l,flow_comm,ierr)
    call mpi_group_free(group_l,ierr)
    if(rank.lt.leader)then
      call mpi_comm_size(flow_comm,ntask,ierr)
    endif
  elseif(leader.eq.0)then
    ! single shared memory partition (nodes=1)
    call mpi_comm_dup(part_comm,flow_comm,ierr)
    call mpi_comm_size(flow_comm,ntask,ierr)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! create communication communicator
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(leader.gt.0)then
    ! at least 2 shared memory partitions (nodes>1)
    allocate(list_rank(leader+1))
    list_rank(1)=0
    do i=2,leader+1
      list_rank(i)=list_rank(i-1)+1
    enddo
    call mpi_group_incl(group_w,flow_comm_lim+1,list_rank,group_l,ierr)
    deallocate(list_rank)
    call mpi_comm_create(mpi_comm_world,group_l,comm_comm,ierr)
    call mpi_group_free(group_l,ierr)
  elseif(leader.eq.0)then
    ! single shared memory partition (nodes=1)
    call mpi_comm_dup(part_comm,comm_comm,ierr)
  endif
#else
  ! if no particles are used, duplicate mpi_comm_world
  call mpi_comm_dup(mpi_comm_world,flow_comm,ierr)
  ntask=ntask_gl
#endif

! rank<flow_comm_lim belongs to flow_comm
if(nodes.gt.1)then
  flow_comm_lim=leader
elseif(nodes.eq.1)then
  flow_comm_lim=ntask_sh
endif


return
end subroutine
