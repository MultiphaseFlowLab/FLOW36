program read_to_paraview

use mpi
use commondata
use wavenumber

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

if(upflag.eq.1)then
  allocate(up(nxf,nzf,nyf))
  allocate(vp(nxf,nzf,nyf))
  allocate(wp(nxf,nzf,nyf))
endif

if(vorflag.eq.1)then
  allocate(omx(nxf,nzf,nyf))
  allocate(omy(nxf,nzf,nyf))
  allocate(omz(nxf,nzf,nyf))
endif

if(strflag.eq.1)then
  allocate(strx(nxf,nzf,nyf))
  allocate(stry(nxf,nzf,nyf))
  allocate(strz(nxf,nzf,nyf))
endif

if(topflag.eq.1)then
  allocate(Qtop(nxf,nzf,nyf))
endif

if(marflag.eq.1)then
  allocate(marx(nxf,nzf,nyf))
  allocate(mary(nxf,nzf,nyf))
  allocate(marz(nxf,nzf,nyf))
endif

if(partposflag.eq.1)then
 allocate(xpar(part_number,3))
 if(partvelflag.eq.1)then
  allocate(upar(part_number,3))
 endif
endif

allocate(kx(nxf/2+1))
allocate(ky(nyf))

kx(1)=0.0d0
do i=2,nxf/2+1
  kx(i)=dble(i-1)*2.0d0*pi/xl
enddo

ky(1)=0.0d0
do i=2,nyf/2+1
  ky(nyf-i+2)=-dble(i-1)*2.0d0*pi/yl
  ky(i)=dble(i-1)*2.0d0*pi/yl
enddo

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

if(upflag.eq.1)then
  deallocate(up,vp,wp)
endif

if(vorflag.eq.1)then
  deallocate(omx,omy,omz)
endif

if(strflag.eq.1)then
  deallocate(strx,stry,strz)
endif

if(topflag.eq.1)then
  deallocate(Qtop)
endif

if(marflag.eq.1)then
  deallocate(marx,mary,marz)
endif

if(partposflag.eq.1)then
 deallocate(xpar)
 if(partvelflag.eq.1) deallocate(upar)
endif

deallocate(kx,ky)

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
