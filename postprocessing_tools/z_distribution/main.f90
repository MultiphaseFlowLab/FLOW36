program mass_center

use mpi
use commondata
use velocity
implicit none

integer :: ierr,i


call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

allocate(phi(nx,nz,ny))
allocate(phic(nx/2+1,nz,ny,2))

call read_grid

call create_plan

do i=nstart,nend,dump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8)') 'Step ',i,' out of ',nend
 call read_fields(i)
enddo


deallocate(x)
deallocate(y)
deallocate(z)
deallocate(phi)
deallocate(phic)

call destroy_plan

call mpi_finalize(ierr)

end program mass_center
