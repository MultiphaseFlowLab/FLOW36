program mass_center

use mpi
use commondata
use wavenumbers
use velocity
use sterm
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
allocate(s1(nx,nz,ny))
allocate(s2(nx,nz,ny))
allocate(s3(nx,nz,ny))

call read_grid

call create_plan

call wavenumber

! create output file
open(42,file='./output/drop_count.dat',status='new',form='formatted')
 write(42,'(3(a16,2x))') 'iteration','t^+','drop count'
close(42,status='keep')

!open(66,file='./output/mass_center.dat',form='formatted',status='new')
! write(66,'(a8,a6,9(a16))') 'nstep','drop','xg','yg','zg','Volume','Area','Mx','My','Mz','(wall units)'
!close(66,status='keep')

!open(67,file='./output/inertia_tensor.dat',form='formatted',status='new')
! write(67,'(2(a6),6(a16))') 'nstep','drop','a11','a12','a13','a22','a23','a33'
!close(67,status='keep')

!open(68,file='./output/eigen_problem.dat',form='formatted',status='new')
! write(68,'(2(a6),4(a16))') 'nstep','drop','eigenvalue','eigenvector x','eigenvector y','eigenvector z'
!close(68,status='keep')

do i=nstart,nend,dump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8)') 'Step ',i,' out of ',nend
 call read_fields(i)
enddo




deallocate(x)
deallocate(y)
deallocate(z)
deallocate(phi)
deallocate(phic)
deallocate(s1)
deallocate(s2)
deallocate(s3)
deallocate(kx)
deallocate(ky)
deallocate(k2)

call destroy_plan

call mpi_finalize(ierr)

end program mass_center
