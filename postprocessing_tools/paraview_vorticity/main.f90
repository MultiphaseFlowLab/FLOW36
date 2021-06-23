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


allocate(u(nx,nz,ny))
allocate(v(nx,nz,ny))
allocate(w(nx,nz,ny))
if(phiflag.eq.1)then
 allocate(phi(nx,nz,ny))
 if(spectral.eq.1)then
  allocate(phic(nx/2+1,nz,ny,2))
 endif
endif
if(psiflag.eq.1)then
 allocate(psi(nx,nz,ny))
 if(spectral.eq.1)then
  allocate(psic(nx/2+1,nz,ny,2))
 endif
endif
if(tempflag.eq.1)then
 allocate(theta(nx,nz,ny))
 if(spectral.eq.1)then
  allocate(thetac(nx/2+1,nz,ny,2))
 endif
endif


!if(spectral.eq.1)then
 allocate(uc(nx/2+1,nz,ny,2))
 allocate(vc(nx/2+1,nz,ny,2))
 allocate(wc(nx/2+1,nz,ny,2))
!endif

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))


call read_grid

!if(spectral.eq.1)then
 call create_plan
!endif

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


deallocate(kx)
deallocate(ky)

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

if(spectral.eq.1)then
 call destroy_plan
endif

call mpi_finalize(ierr)

end program read_to_paraview
