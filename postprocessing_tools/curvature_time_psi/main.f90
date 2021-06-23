program curvature

use mpi
use commondata
! use wavenumbers
use pdf_calc
implicit none

integer :: ierr,i
double precision :: dd


call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(phi(nx,nz,ny))
allocate(kv(nx,nz,ny))
allocate(phic(nx/2+1,nz,ny,2))
allocate(psi(nx*exp_x,(nz-1)*exp_z+1,ny*exp_y))

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
call read_grid

allocate(xfg(nx*exp_x))
allocate(yfg(ny*exp_y))
allocate(zfg((nz-1)*exp_z+1))
xfg(1)=0.0d0
dd=xl/dble(nx*exp_x-1)
do i=2,exp_x*nx
  xfg(i)=xfg(i-1)+dd
enddo
yfg(1)=0.0d0
dd=yl/dble(ny*exp_y-1)
do i=2,exp_y*ny
  yfg(i)=yfg(i-1)+dd
enddo
do i=1,(nz-1)*exp_z+1
  zfg(i)=dcos((dble(i-1)*pi)/(((nz-1)*exp_z+1)-1))
enddo



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
deallocate(psi)
deallocate(kv)
deallocate(x,y,z,xfg,yfg,zfg)

call destroy_plan

call mpi_finalize(ierr)

end program curvature
