program main

use mpi
use commondata
use fields
use wavenumber

integer :: i,j

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(u(nx,nz,ny))
allocate(v(nx,nz,ny))
allocate(w(nx,nz,ny))
allocate(phi(nx,nz,ny))
if(phi_flag.eq.0) phi=-1.0d0
if(psi_flag.eq.1) then
 allocate(psi(expx*nx,expz*(nz-1)+1,expy*ny))
endif

allocate(kx(nx/2+1))
allocate(ky(ny))
allocate(k2(nx/2+1,ny))
allocate(kxfg(nxfg/2+1))
allocate(kyfg(nyfg))
allocate(k2fg(nxfg/2+1,nyfg))

! wavenumbers
kx(1)=0.0d0
do i=2,nx/2+1
  kx(i)=dble(i-1)*2.0d0*pi/xl
enddo

ky(1)=0.0d0
do i=2,ny/2+1
  ky(ny-i+2)=-dble(i-1)*2.0d0*pi/yl
  ky(i)=dble(i-1)*2.0d0*pi/yl
enddo

do j=1,ny
  do i=1,nx/2+1
    k2(i,j)=kx(i)*kx(i)+ky(j)*ky(j)
  enddo
enddo

kxfg(1)=0.0d0
do i=2,nxfg/2+1
  kxfg(i)=dble(i-1)*2.0d0*pi/xl
enddo

kyfg(1)=0.0d0
do i=2,nyfg/2+1
  kyfg(nyfg-i+2)=-dble(i-1)*2.0d0*pi/yl
  kyfg(i)=dble(i-1)*2.0d0*pi/yl
enddo

do j=1,nyfg
  do i=1,nxfg/2+1
    k2fg(i,j)=kxfg(i)*kxfg(i)+kyfg(j)*kyfg(j)
  enddo
enddo

! plan creation
call create_plan
call create_plan_fg
call create_plan_2D

call read_grid

do i=nstart,nend,delta
 call marangoni(i)
enddo
! run last step if not already run
if(mod(nend,delta).ne.0) call marangoni(nend)



call destroy_plan

deallocate(x,y,z)
deallocate(u,v,w)
deallocate(phi)
if(psi_flag.eq.1) deallocate(psi)
deallocate(kx,ky,k2,kxfg,kyfg,k2fg)


call mpi_finalize(ierr)

return
end program main
