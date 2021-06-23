program psi_V

use mpi
use commondata
use wavenumbers
use sterm
implicit none

integer :: ierr,i
double precision :: dx,dy


call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(xfg(nx*exp_x))
allocate(yfg(ny*exp_y))
allocate(zfg((nz-1)*exp_z+1))

! define xfg,yfg,zfg on fine grid
dx=xl/dble(exp_x*nx-1)
dy=yl/dble(exp_y*ny-1)

do i=1,exp_x*nx
  xfg(i)=dble(i-1)*dx
enddo

do i=1,exp_y*ny
  yfg(i)=dble(i-1)*dy
enddo

do i=1,exp_z*(nz-1)+1
  zfg(i)=dcos(((i-1)*pi)/((exp_z*(nz-1)+1)-1))
enddo


allocate(phi(nx,nz,ny))
allocate(psi(nx*exp_x,(nz-1)*exp_z+1,ny*exp_y))


call read_grid

call create_plan
call wavenumber

do i=nstart,nend,dump
 if(rank.eq.0) write(*,'(1x,a,i8,a,i8)') 'Step ',i,' out of ',nend
 call read_fields(i)
enddo




deallocate(x,y,z)
deallocate(xfg,yfg,zfg)
deallocate(phi)
deallocate(psi)
deallocate(kx)
deallocate(ky)
deallocate(k2)

call destroy_plan

call mpi_finalize(ierr)

end program psi_V
