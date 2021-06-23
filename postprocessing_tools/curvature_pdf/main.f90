program curvature

use mpi
use commondata
use wavenumbers
use pdf_calc
implicit none

integer :: ierr,i


call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(phi(nx,nz,ny))
allocate(kv(nx,nz,ny))
allocate(phic(nx/2+1,nz,ny,2))


allocate(x(nx))
allocate(y(ny))
allocate(z(nz)) ! deallocate missing
call read_grid

call create_plan

call wavenumber

mink=0.0d0
maxk=0.0d0

do i=nstart,nend,ndump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8,a)') 'Step ',i,' out of ',nend,' , calculating curvature'
 call read_fields(i)
enddo

! apply regola della serranda
mink=mink-abs(mink)*0.1d0
maxk=maxk+abs(maxk)*0.1d0

do i=nstart,nend,ndump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8,a)') 'Step ',i,' out of ',nend,' , generating curvature pdf'
 call read_curvature(i)
enddo



deallocate(phi)
deallocate(phic)
deallocate(kv)
deallocate(kx)
deallocate(ky)
deallocate(k2)

call destroy_plan

call mpi_finalize(ierr)

end program curvature
