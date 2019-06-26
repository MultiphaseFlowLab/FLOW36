subroutine split_comm

use mpi
use commondata

#define particles particlecompflag


#if particles == 1
! split flow, particle and communication communicators
! assign rank (global), ntask_gl (global) and ntask (local)

! duplicate mpi_comm_world, only temporary until particle part done
call mpi_comm_dup(mpi_comm_world,flow_comm,ierr)

! check that there are more than 1 nodes in case of particles
! eventually flow_comm and part_comm might be the same


#else
  ! if no particles are used, duplicate mpi_comm_world
  call mpi_comm_dup(mpi_comm_world,flow_comm,ierr)
#endif

! query number of MPI processes, assign rank number
call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

return
end subroutine
