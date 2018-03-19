program statistics

use mpi
use commondata
implicit none

double precision, allocatable :: stats(:,:),bufr(:)

integer :: ierr,i,dump,nstep,k,iadd

character(len=5) :: namevar

call mpi_init(ierr)

call mpi_comm_size(mpi_comm_world,ntask,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)

call read_input
call mpi_barrier(mpi_comm_world,ierr)

allocate(u(nx,nz,ny))
allocate(mean_u(nz))
allocate(rms_u(nz))
allocate(skw_u(nz))
allocate(flt_u(nz))

mean_u=0.0d0
rms_u=0.0d0
skw_u=0.0d0
flt_u=0.0d0


if(spectral.eq.1)then
 allocate(uc(nx/2+1,nz,ny,2))
endif

allocate(x(nx))
allocate(y(ny))
allocate(z(nz))


call read_grid

call create_plan

if(spectral.eq.1)then
 dump=sdump
else
 dump=ndump
endif

! rank 0: u stats
! rank 1: v stats
! rank 3: w stats
if(rank.eq.0)then
 if(spectral.eq.1)then
  namevar='uc_  '
 else
  namevar='u_   '
 endif
elseif(rank.eq.1)then
 if(spectral.eq.1)then
  namevar='vc_  '
 else
  namevar='v_   '
 endif
elseif(rank.eq.2)then
 if(spectral.eq.1)then
  namevar='wc_  '
 else
  namevar='w_   '
 endif
endif

counter=0
nstep=0

! start mean velocity calculation
do i=nstart,nend,dump
 nstep=i
 call read_fields(i,namevar)
enddo
! open file nend, unless already opened in previous loop
if(nstep.ne.nend) then
 call read_fields(nend,namevar)
endif

mean_u=mean_u/dble(nx*ny*counter)

rms_u=rms_u/dble(nx*ny*counter)
rms_u=rms_u**0.5d0

skw_u=skw_u/dble(nx*ny*counter)
flt_u=flt_u/dble(nx*ny*counter)

do k=2,nz-1
  skw_u(k)=skw_u(k)/rms_u(k)**3
  flt_u(k)=flt_u(k)/rms_u(k)**4
enddo
skw_u(1)=0.0d0
flt_u(1)=0.0d0
skw_u(nz)=0.0d0
flt_u(nz)=0.0d0


if(rank.eq.0)write(*,*) 'Gathering results to rank 0'
! gather results
if(rank.ne.0)then
! send to 0
! mean
call mpi_send(mean_u,nz,mpi_double_precision,0,rank*10+1,mpi_comm_world,ierr)
! rms
call mpi_send(rms_u,nz,mpi_double_precision,0,rank*10+2,mpi_comm_world,ierr)
! skw
call mpi_send(skw_u,nz,mpi_double_precision,0,rank*10+3,mpi_comm_world,ierr)
! flt
call mpi_send(flt_u,nz,mpi_double_precision,0,rank*10+4,mpi_comm_world,ierr)

else
! z, u_m, v_m, w_m, rms_u, rms_v, rms_w, skw_u, skw_v, skw_w, flt_u, flt_v, flt_w
allocate(stats(nz,13))
allocate(bufr(nz))
stats(:,1)=z(:)
stats(:,2)=mean_u(:)
stats(:,5)=rms_u(:)
stats(:,8)=skw_u(:)
stats(:,11)=flt_u(:)
! receive from other ranks
! receive v mean
call mpi_recv(bufr,nz,mpi_double_precision,1,11,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,3)=bufr(:)
! receive w mean
call mpi_recv(bufr,nz,mpi_double_precision,2,21,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,4)=bufr(:)
! receive v rms
call mpi_recv(bufr,nz,mpi_double_precision,1,12,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,6)=bufr(:)
! receive w rms
call mpi_recv(bufr,nz,mpi_double_precision,2,22,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,7)=bufr(:)
! receive v skw
call mpi_recv(bufr,nz,mpi_double_precision,1,13,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,9)=bufr(:)
! receive w skw
call mpi_recv(bufr,nz,mpi_double_precision,2,23,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,10)=bufr(:)
! receive v flt
call mpi_recv(bufr,nz,mpi_double_precision,1,14,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,12)=bufr(:)
! receive w flt
call mpi_recv(bufr,nz,mpi_double_precision,2,24,mpi_comm_world,mpi_status_ignore,ierr)
if(.not.mpi_async_protects_nonblocking) call mpi_get_address(bufr,iadd,ierr)
stats(:,13)=bufr(:)

call write_output(stats)

deallocate(bufr)
deallocate(stats)
endif




deallocate(u)
deallocate(mean_u)
deallocate(rms_u)
deallocate(skw_u)
deallocate(flt_u)

if(spectral.eq.1)then
 deallocate(uc)
endif

deallocate(x)
deallocate(y)
deallocate(z)

call destroy_plan

call mpi_finalize(ierr)

end program statistics
