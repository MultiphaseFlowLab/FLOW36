subroutine read_fields(u,nt,namevar,restart)

use commondata
use mpi
use mpiIO
use par_size

double precision :: u(nx,fpz,fpy)

integer :: nt,restart
integer :: f_handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
logical :: check
integer :: file_size

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
endif

offset=0

inquire(file=trim(fname),exist=check,size=file_size)
if(file_size.ne.8*nx*ny*nz)then
  if(rank.eq.0) write(*,*) 'Wrong file size'
  call exit(0)
endif

if(check.eqv..true.)then
 call mpi_file_open(flow_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'internal',mpi_info_null,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)
 call mpi_file_read(f_handle,u,nx*fpy*fpz,mpi_double_precision,mpi_status_ignore,ierr)
 call mpi_file_close(f_handle,ierr)
else
 ! stop simulation
 if(rank.eq.0) write(*,*) 'Missing ',trim(namevar),' fields'
 call exit(0)
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_fields_s(u,nt,namevar,restart)

use commondata
use mpi
use mpiIO
use par_size
use shrink_grid

double precision :: u(spx,nz,spy,2),ucd(dimc(1),dimc(2),dimc(3),2)

integer :: nt
integer :: f_handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
logical :: check
integer :: file_size,mx,my,mz

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
endif

offset=0

mx=floor(2.0d0/3.0d0*dble(nx/2+1))
my=floor(2.0d0/3.0d0*dble(ny/2+1))+floor(2.0d0/3.0d0*dble(ny/2))
mz=floor(2.0d0/3.0d0*dble(nz))
inquire(file=trim(fname),exist=check,size=file_size)
if(file_size.ne.8*mx*my*mz*2)then
  if(rank.eq.0) write(*,*) 'Wrong file size (spectral)'
  call exit(0)
endif

if(check.eqv..true.)then
 call mpi_file_open(sp_save_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'external32',mpi_info_null,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'internal',mpi_info_null,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'native',mpi_info_null,ierr)
 call mpi_file_read(f_handle,ucd,dimc(1)*dimc(2)*dimc(3)*2,mpi_double_precision,mpi_status_ignore,ierr)
 call mpi_file_close(f_handle,ierr)
 call expand_domain(ucd,u)
else
 ! stop simulation
 if(rank.eq.0) write(*,*) 'Missing ',trim(namevar),' fields'
 call exit(0)
endif

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_fields_fg(u,nt,namevar,restart)

use commondata
use mpi
use mpiIO
use par_size
use dual_grid

double precision :: u(npsix,fpzpsi,fpypsi)

integer :: nt,restart
integer :: f_handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
logical :: check
integer :: file_size

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'_fg.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_fg_'//time//'.dat'
endif

offset=0

inquire(file=trim(fname),exist=check,size=file_size)
if(file_size.ne.8*npsix*npsiy*npsiz)then
  if(rank.eq.0) write(*,*) 'Wrong file size in fine grid'
  call exit(0)
endif

if(check.eqv..true.)then
 call mpi_file_open(flow_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'external32',mpi_info_null,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'internal',mpi_info_null,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'native',mpi_info_null,ierr)
 call mpi_file_read(f_handle,u,npsix*fpypsi*fpzpsi,mpi_double_precision,mpi_status_ignore,ierr)
 call mpi_file_close(f_handle,ierr)
else
 ! stop simulation
 if(rank.eq.0) write(*,*) 'Missing ',trim(namevar),' fields'
 call exit(0)
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_fields_s_fg(u,nt,namevar,restart)

use commondata
use mpi
use mpiIO
use par_size
use dual_grid
use shrink_grid

double precision :: u(spxpsi,npsiz,spypsi,2),ucd(dimc_fg(1),dimc_fg(2),dimc_fg(3),2)

integer :: nt
integer :: f_handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
logical :: check
integer :: file_size,mx,my,mz

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'_fg.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_fg_'//time//'.dat'
endif

offset=0

mx=floor(2.0d0/3.0d0*dble(npsix/2+1))
my=floor(2.0d0/3.0d0*dble(npsiy/2+1))+floor(2.0d0/3.0d0*dble(npsiy/2))
mz=floor(2.0d0/3.0d0*dble(npsiz))
inquire(file=trim(fname),exist=check,size=file_size)
if(file_size.ne.8*mx*my*mz*2)then
  if(rank.eq.0) write(*,*) 'Wrong file size in fine grid (spectral)'
  call exit(0)
endif

if(check.eqv..true.)then
 call mpi_file_open(sp_save_comm_fg,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype_fg,'native',mpi_info_null,ierr)
 call mpi_file_read(f_handle,ucd,dimc_fg(1)*dimc_fg(2)*dimc_fg(3)*2,mpi_double_precision,mpi_status_ignore,ierr)
 call mpi_file_close(f_handle,ierr)
 call expand_domain_fg(ucd,u)
else
 ! stop simulation
 if(rank.eq.0) write(*,*) 'Missing ',trim(namevar),' fields'
 call exit(0)
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_fields_serial(u,nt,namevar,restart)

use commondata
use mpi
use mpiIO
use par_size

double precision :: u(nx,fpz,fpy),temp_u(nx,nz,ny)

integer :: nt
integer :: file_size

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar

logical :: check

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
endif

inquire(file=trim(fname),exist=check,size=file_size)
if(file_size.ne.8*nx*ny*nz)then
  if(rank.eq.0) write(*,*) 'Wrong file size'
  call exit(0)
endif

if(check.eqv..true.)then
 open(unit=66,file=fname,form='unformatted',status='old',access='stream')
 read(66) temp_u
 close(66,status='keep')
else
 ! stop simulation
 if(rank.eq.0) write(*,*) 'Missing ',trim(namevar),' fields'
 call exit(0)
endif

u(1:nx,1:fpz,1:fpy)=temp_u(1:nx,fstart(2)+1:fstart(2)+fpz,fstart(3)+1:fstart(3)+fpy)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_part_pos(nt,restart)

use commondata
use mpi
use mpiIO
use particle

double precision :: varloc(part_index(rank_loc+1,2),3)

integer :: nt,restart,j
integer :: f_handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=3) :: setnum
logical :: check

do j=1,nset
  write(setnum,'(i3.3)') j

  if(restart.eq.0)then
   fname='./initial_fields/pos_'//setnum//'.dat'
  else
   write(time,'(I8.8)') nt
   fname=trim(folder)//'/pos_'//setnum//'_'//time//'.dat'
  endif

  offset=0

  inquire(file=trim(fname),exist=check)

  if(check.eqv..true.)then
   call mpi_file_open(part_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
   !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
   !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
   call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)
   call mpi_file_read(f_handle,varloc,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)
   call mpi_file_close(f_handle,ierr)
   xp(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:,j)=varloc
  else
   ! stop simulation
   if(rank.eq.leader) write(*,*) 'Missing pospar fields'
   call exit(0)
  endif
enddo

call mpi_win_fence(0,window_xp,ierr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_part_vel(nt,restart)

use commondata
use mpi
use mpiIO
use particle

double precision :: varloc(part_index(rank_loc+1,2),3)

integer :: nt,restart,j
integer :: f_handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=3) :: setnum
logical :: check

do j=1,nset
  write(setnum,'(i3.3)') j

  if(restart.eq.0)then
   fname='./initial_fields/vel_'//setnum//'.dat'
  else
   write(time,'(I8.8)') nt
   fname=trim(folder)//'/vel_'//setnum//'_'//time//'.dat'
  endif

  offset=0

  inquire(file=trim(fname),exist=check)

  if(check.eqv..true.)then
   call mpi_file_open(part_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)
   !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
   !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
   call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)
   call mpi_file_read(f_handle,varloc,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)
   call mpi_file_close(f_handle,ierr)
   up(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:,j)=varloc
  else
   ! stop simulation
   if(rank.eq.leader) write(*,*) 'Missing velpar fields'
   call exit(0)
  endif
enddo

call mpi_win_fence(0,window_up,ierr)

return
end
