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

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
endif

offset=0

call mpi_file_open(flow_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'internal',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

call mpi_file_read_all(f_handle,u,nx*fpy*fpz,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

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

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
endif

offset=0

call mpi_file_open(sp_save_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'internal',mpi_info_null,ierr)

call mpi_file_read_all(f_handle,ucd,dimc(1)*dimc(2)*dimc(3)*2,mpi_double_precision,mpi_status_ignore,ierr)

call mpi_file_close(f_handle,ierr)

call expand_domain(ucd,u)

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

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'_fg.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_fg_'//time//'.dat'
endif

offset=0

call mpi_file_open(flow_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'external32',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'internal',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'native',mpi_info_null,ierr)

call mpi_file_read_all(f_handle,u,npsix*fpypsi*fpzpsi,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

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

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'_fg.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_fg_'//time//'.dat'
endif

offset=0

call mpi_file_open(sp_save_comm_fg,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype_fg,'internal',mpi_info_null,ierr)

call mpi_file_read_all(f_handle,ucd,dimc_fg(1)*dimc_fg(2)*dimc_fg(3)*2,mpi_double_precision,mpi_status_ignore,ierr)

call mpi_file_close(f_handle,ierr)

call expand_domain_fg(ucd,u)

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


character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar

if(restart.eq.0)then
 fname='./initial_fields/'//trim(namevar)//'.dat'
else
 write(time,'(I8.8)') nt
 fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
endif


 open(unit=66,file=fname,form='unformatted',status='old',access='stream')

 read(66) temp_u

 close(66,status='keep')


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

do j=1,nset
  write(setnum,'(i3.3)') j

  if(restart.eq.0)then
   fname='./initial_fields/pos_'//setnum//'.dat'
  else
   write(time,'(I8.8)') nt
   fname=trim(folder)//'/pos_'//setnum//'_'//time//'.dat'
  endif

  offset=0

  call mpi_file_open(part_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

  !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
  !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)

  call mpi_file_read_all(f_handle,varloc,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)

  call mpi_file_close(f_handle,ierr)

  xp(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:,j)=varloc
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


do j=1,nset
  write(setnum,'(i3.3)') j

  if(restart.eq.0)then
   fname='./initial_fields/vel_'//setnum//'.dat'
  else
   write(time,'(I8.8)') nt
   fname=trim(folder)//'/vel_'//setnum//'_'//time//'.dat'
  endif

  offset=0

  call mpi_file_open(part_comm,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

  !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
  !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)

  call mpi_file_read_all(f_handle,varloc,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)

  call mpi_file_close(f_handle,ierr)

  up(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:,j)=varloc
enddo

call mpi_win_fence(0,window_up,ierr)

return
end
