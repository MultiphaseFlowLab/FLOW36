program invariants

use mpi
use commondata
use wavenumber

implicit none

integer :: ierr,i,dump

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(u(nx,nz,ny))
allocate(v(nx,nz,ny))
allocate(w(nx,nz,ny))
allocate(uc(nx/2+1,nz,ny,2))
allocate(vc(nx/2+1,nz,ny,2))
allocate(wc(nx/2+1,nz,ny,2))

allocate(phi(nx,nz,ny))


allocate(kx(nx/2+1))
allocate(ky(ny))

kx(1)=0.0d0
do i=2,nx/2+1
  kx(i)=dble(i-1)*2.0d0*pi/xl
enddo

ky(1)=0.0d0
do i=2,ny/2+1
  ky(ny-i+2)=-dble(i-1)*2.0d0*pi/yl
  ky(i)=dble(i-1)*2.0d0*pi/yl
enddo

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

call read_grid

call create_plan
call create_plan_fg

if(spectral.eq.1)then
 dump=sdump
else
 dump=ndump
endif

! samples_in=0
! samples_out=0

call initialize_jpdf

do i=nstart,nend,dump
 call read_fields(i)
enddo

! write(*,*) samples_in,samples_out
! call join_samples

call write_jpdf

call clear_jpdf


deallocate(u,v,w)
deallocate(uc,vc,wc)

deallocate(phi)

deallocate(kx,ky)

deallocate(x)
deallocate(y)
deallocate(z)

call destroy_plan

call mpi_finalize(ierr)

end program invariants
