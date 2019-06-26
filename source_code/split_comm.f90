subroutine split_comm

use mpi
use commondata

#define particles particlescompflag

! split flow, particle and communication communicators

! assign rank (global), ntask and ntask local

! duplicate mpi_comm_world, only temporary until particle part done
! or in case no particles are used
call mpi_comm_dup(mpi_comm_world,flow_comm,ierr)


! query number of MPI processes, assign rank number
call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

return
end subroutine
