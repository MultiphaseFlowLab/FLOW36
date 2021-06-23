subroutine read_fields(u,nt,namevar)

use commondata
use mpi
use mpiIO

double precision :: u(nx,fpz,fpy)

integer :: nt
integer :: f_handle

integer(mpi_offset_kind) :: offset

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar

write(time,'(I8.8)') nt
fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'

offset=0

call mpi_file_open(mpi_comm_world,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'internal',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

call mpi_file_read_all(f_handle,u,nx*fpy*fpz,mpi_double_precision,mpi_status_ignore,ierr)

call mpi_file_close(f_handle,ierr)

return
end





subroutine read_fields_s(u,nt,namevar)

use commondata
use mpi
use mpiIO

double precision :: u(spx,nz,spy,2)

integer :: nt
integer :: f_handle

integer(mpi_offset_kind) :: offset

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar


write(time,'(I8.8)') nt
fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'

offset=0

call mpi_file_open(mpi_comm_world,fname,mpi_mode_rdonly,mpi_info_null,f_handle,ierr)

call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'internal',mpi_info_null,ierr)

call mpi_file_read_all(f_handle,u,spx*spy*nz*2,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

return
end



! 
!
! subroutine read_fields_serial(u,nt,namevar,restart)
!
! use commondata
! use mpi
! use mpiIO
!
! double precision :: u(nx,fpz,fpy),temp_u(nx,nz,ny)
!
! integer :: nt
!
!
! character(len=8) :: time
! character(len=40) :: fname
! character(len=5) :: namevar
!
! if(restart.eq.0)then
!  fname='./initial_fields/'//trim(namevar)//'.dat'
! else
!  write(time,'(I8.8)') nt
!  fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'
! endif
!
!
!  open(unit=66,file=fname,form='unformatted',status='old',access='stream')
!
!  read(66) temp_u
!
!  close(66,status='keep')
!
!
! u(1:nx,1:fpz,1:fpy)=temp_u(1:nx,fstart(2)+1:fstart(2)+fpz,fstart(3)+1:fstart(3)+fpy)
!
!
! return
! end
