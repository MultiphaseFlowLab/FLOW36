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

allocate(u(nxf,nzf,nyf))
allocate(v(nxf,nzf,nyf))
allocate(w(nxf,nzf,nyf))
allocate(uc(nxf/2+1,nzf,nyf,2))
allocate(vc(nxf/2+1,nzf,nyf,2))
allocate(wc(nxf/2+1,nzf,nyf,2))

if(phiflag.eq.1)then
 allocate(phi(nxf,nzf,nyf))
 allocate(phic(nxf/2+1,nzf,nyf,2))
endif
if(psiflag.eq.1)then
 allocate(psi(nxf,nzf,nyf))
 allocate(psic(nxf/2+1,nzf,nyf,2))
endif
if(tempflag.eq.1)then
 allocate(theta(nxf,nzf,nyf))
 allocate(thetac(nxf/2+1,nzf,nyf,2))
endif



allocate(x(nx))
allocate(y(ny))
allocate(z(nz))

allocate(xfg(nxf))
allocate(yfg(nyf))
allocate(zfg(nzf))


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

if(spectral.eq.1)then
 deallocate(uc)
 deallocate(vc)
 deallocate(wc)
endif

if(phiflag.eq.1)then
 deallocate(phi)
 if(spectral.eq.1)then
  deallocate(phic)
 endif
endif

if(psiflag.eq.1)then
 deallocate(psi)
 if(spectral.eq.1)then
  deallocate(psic)
 endif
endif

if(tempflag.eq.1)then
 deallocate(theta)
 if(spectral.eq.1)then
  deallocate(thetac)
 endif
endif

deallocate(x)
deallocate(y)
deallocate(z)

deallocate(xfg)
deallocate(yfg)
deallocate(zfg)


if(spectral.eq.1)then
 call destroy_plan
endif

call mpi_finalize(ierr)

end program read_to_paraview