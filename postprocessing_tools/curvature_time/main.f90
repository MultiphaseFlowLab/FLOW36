program curvature

use mpi
use commondata
! use wavenumbers
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
allocate(z(nz))
call read_grid

call create_plan
! call wavenumber

fhandle=666
open(fhandle,file='./output/curvature.dat',status='new',form='formatted')
write(fhandle,'(4(a16,2x))') 'step','t^+','mean k','rms k'

do i=nstart,nend,ndump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8,a)') 'Step ',i,' out of ',nend,' , calculating curvature'
 call read_fields(i)
enddo

close(fhandle,status='keep')

deallocate(phi)
deallocate(phic)
deallocate(kv)
deallocate(x,y,z)

call destroy_plan

call mpi_finalize(ierr)

end program curvature
