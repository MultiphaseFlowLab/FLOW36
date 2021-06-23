program read_to_paraview

use mpi
use commondata
implicit none

integer :: ierr,i,nstep,dump

logical :: check

character(len=40) :: namefile

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(u(nx,nz,ny))
allocate(v(nx,nz,ny))
allocate(w(nx,nz,ny))

allocate(u_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny))
allocate(v_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny))
allocate(w_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny))

if(phiflag.eq.1)then
 allocate(phi(nx,nz,ny))
 allocate(phi_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny))
endif



if(psiflag.eq.1)then
 allocate(psi_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny))
endif
if(tempflag.eq.1)then
 allocate(theta(nx,nz,ny))
 allocate(theta_fg(exp_x*nx,exp_z*(nz-1)+1,exp_y*ny))
endif


allocate(x(nx))
allocate(y(ny))
allocate(z(nz))
allocate(x_fg(nx*exp_x))
allocate(y_fg(ny*exp_y))
allocate(z_fg((nz-1)*exp_z+1))


call read_grid

call create_plan
call create_plan_fg

if(spectral.eq.1)then
 dump=sdump
else
 dump=ndump
endif

do i=nstart,nend,ntask*dump
 nstep=i+rank*dump
! write(*,*) rank,i,nstep
 call read_fields(nstep)
enddo

! barrier to be sure that there isn't another rank handling that file but
! it still hasn't finished yet
call mpi_barrier(mpi_comm_world,ierr)
! open file nend
if(rank.eq.0)then
 write(namefile,'(a,i8.8,a)') './output/OUTPAR_',nend,'.vtk'
 inquire(file=trim(namefile),exist=check)
 if(check.eqv..false.)then
  call read_fields(nend)
 endif
endif



deallocate(u)
deallocate(v)
deallocate(w)
deallocate(u_fg,v_fg,w_fg)


if(phiflag.eq.1)then
 deallocate(phi)
endif

if(psiflag.eq.1)then
 deallocate(psi_fg)
endif

if(tempflag.eq.1)then
 deallocate(theta)
 deallocate(theta_fg)
endif

deallocate(x)
deallocate(y)
deallocate(z)
deallocate(x_fg,y_fg,z_fg)

if(spectral.eq.1)then
 call destroy_plan
endif

call mpi_finalize(ierr)

end program read_to_paraview
