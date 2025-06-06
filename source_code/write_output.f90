subroutine write_output(u,nt,namevar)

use mpi
use commondata
use par_size
use mpiIo

double precision :: u(nx,fpz,fpy)

integer :: nt
integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar

write(time,'(I8.8)') nt

fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'

offset=0

call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

call mpi_file_write_all(f_handle,u,nx*fpy*fpz,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_spectral(uc,nt,namevar)

use mpi
use commondata
use par_size
use mpiIo
use shrink_grid

double precision :: uc(spx,nz,spy,2),ucd(dimc(1),dimc(2),dimc(3),2)

integer :: nt
integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar


call shrink_domain(uc,ucd)

write(time,'(I8.8)') nt

fname=trim(folder)//'/'//trim(namevar)//'_'//time//'.dat'

offset=0

! call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
call mpi_file_open(sp_save_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

! call mpi_file_write_all(f_handle,uc,spx*spy*nz*2,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_write_all(f_handle,ucd,dimc(1)*dimc(2)*dimc(3)*2,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_fg(u,nt,namevar)

use mpi
use commondata
use par_size
use mpiIo
use dual_grid

double precision :: u(npsix,fpzpsi,fpypsi)

integer :: nt
integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar

write(time,'(I8.8)') nt

fname=trim(folder)//'/'//trim(namevar)//'_fg_'//time//'.dat'

offset=0

call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'native',mpi_info_null,ierr)

call mpi_file_write_all(f_handle,u,npsix*fpypsi*fpzpsi,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_spectral_fg(uc,nt,namevar)

use mpi
use commondata
use par_size
use mpiIo
use dual_grid
use shrink_grid

double precision :: uc(spxpsi,npsiz,spypsi,2),ucd(dimc_fg(1),dimc_fg(2),dimc_fg(3),2)

integer :: nt
integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar


call shrink_domain_fg(uc,ucd)

write(time,'(I8.8)') nt

fname=trim(folder)//'/'//trim(namevar)//'_fg_'//time//'.dat'

offset=0

! call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
call mpi_file_open(sp_save_comm_fg,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype_fg,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

! call mpi_file_write_all(f_handle,uc,spxpsi*spypsi*npsiz*2,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_write_all(f_handle,ucd,dimc_fg(1)*dimc_fg(2)*dimc_fg(3)*2,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine write_output_many(u,nt,namevar)

!use commondata
!use par_size

!double precision :: u(nx,fpz,fpy)

!integer :: nt

!character(len=8) :: time, srank
!character(len=40) :: fname
!character(len=5) :: namevar

!write(time,'(I8.8)') nt
!write(srank,'(I8.8)') rank

!fname=trim(folder)//'/'//trim(namevar)//'_'//time//'_r'//srank//'.dat'

! open(66+rank,file=fname,form='unformatted',status='new',convert='little_endian')

! write(66+rank) u

! close(66+rank,status='keep')

!return
!end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_part(var,nt,namevar)

use mpi
use mpiIO
use commondata
use particle

integer :: nt,i
integer :: f_handle ! file handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
character(len=3) :: setnum

double precision :: var(part_number,3,nset),varloc(part_index(rank_loc+1,2),3)

write(time,'(I8.8)') nt

do i=1,nset
 write(setnum,'(i3.3)') i
 fname=trim(folder)//'/'//trim(namevar)//'_'//trim(setnum)//'_'//time//'.dat'

 offset=0
 varloc=var(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:,i)

 call mpi_file_open(part_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)

 call mpi_file_write_all(f_handle,varloc,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)

 call mpi_file_close(f_handle,ierr)
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_partf(nt)

use mpi
use mpiIO
use commondata
use particle

integer :: nt,i,j
integer :: f_handle ! file handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
character(len=3) :: setnum

double precision :: var(part_index(rank_loc+1,2),3)

write(time,'(I8.8)') nt
namevar='vel_f'

do j=1,nset
 write(setnum,'(i3.3)') j
 fname=trim(folder)//'/'//trim(namevar)//'_'//trim(setnum)//'_'//time//'.dat'

 do i=1,part_index(rank_loc+1,2)
   call lagran4(xp(part_index(rank_loc+1,1)+i,:,j),var(i,:))
 enddo

 offset=0

 call mpi_file_open(part_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)

 call mpi_file_write_all(f_handle,var,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)

 call mpi_file_close(f_handle,ierr)
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_partT(nt)

use mpi
use mpiIO
use commondata
use particle

integer :: nt,i,j
integer :: f_handle ! file handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
character(len=3) :: setnum

double precision :: var(part_index(rank_loc+1,2))

write(time,'(I8.8)') nt
namevar='Tp_f'

do j=1,nset
 write(setnum,'(i3.3)') j
 fname=trim(folder)//'/'//trim(namevar)//'_'//trim(setnum)//'_'//time//'.dat'

 do i=1,part_index(rank_loc+1,2)
   call lagran4_T(xp(part_index(rank_loc+1,1)+i,:,j),var(i))
 enddo

 offset=0

 call mpi_file_open(part_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save_scalar,'external32',mpi_info_null,ierr)
 !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save_scalar,'internal',mpi_info_null,ierr)
 call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save_scalar,'native',mpi_info_null,ierr)

 call mpi_file_write_all(f_handle,var,part_index(rank_loc+1,2),mpi_double_precision,mpi_status_ignore,ierr)

 call mpi_file_close(f_handle,ierr)
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_failure(nt)

use commondata
use velocity
use phase_field
use sim_par
use surfactant
use temperature
use particle

use, intrinsic :: ISO_C_BINDING

character(len=80) :: in,out
character(len=3) :: setnum

integer :: nt,ren,j

logical :: check

#define machine machineflag
#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define particles particlecompflag
#define stat_dump stats_dump_frequency
#define spectral_flag savespectralflag

#if machine == 4
interface
  function rename (old_path, new_path) bind (c)
    use iso_c_binding
    integer (kind=c_int) :: rename
    character (kind=c_char) :: old_path(*)
    character (kind=c_char) :: new_path(*)
  end function rename
end interface
#endif


if(rank.eq.0)then
  inquire(file=trim(folder)//'/backup/time_check.dat',exist=check)
#if machine == 4
  if(check.eqv..true.)then
    in=trim(folder)//'/backup/time_check.dat'
    out=trim(folder)//'/backup/time_check_old.dat'
    call copy_append_old(in,out,nt)
    ren=rename(trim(folder)//'/backup/iteration.dat'//C_NULL_CHAR,trim(folder)//'/backup/iteration_old.dat'//C_NULL_CHAR)
#if spectral_flag == 1
    ren=rename(trim(folder)//'/backup/uc.dat'//C_NULL_CHAR,trim(folder)//'/backup/uc_old.dat'//C_NULL_CHAR)
    ren=rename(trim(folder)//'/backup/vc.dat'//C_NULL_CHAR,trim(folder)//'/backup/vc_old.dat'//C_NULL_CHAR)
    ren=rename(trim(folder)//'/backup/wc.dat'//C_NULL_CHAR,trim(folder)//'/backup/wc_old.dat'//C_NULL_CHAR)
#if phiflag == 1
    ren=rename(trim(folder)//'/backup/phic.dat'//C_NULL_CHAR,trim(folder)//'/backup/phic_old.dat'//C_NULL_CHAR)
#if psiflag == 1
    ren=rename(trim(folder)//'/backup/psic_fg.dat'//C_NULL_CHAR,trim(folder)//'/backup/psic_fg_old.dat'//C_NULL_CHAR)
#endif
#endif
#if tempflag == 1
    ren=rename(trim(folder)//'/backup/Tc.dat'//C_NULL_CHAR,trim(folder)//'/backup/Tc_old.dat'//C_NULL_CHAR)
#endif
#elif spectral_flag == 0
    ren=rename(trim(folder)//'/backup/u.dat'//C_NULL_CHAR,trim(folder)//'/backup/u_old.dat'//C_NULL_CHAR)
    ren=rename(trim(folder)//'/backup/v.dat'//C_NULL_CHAR,trim(folder)//'/backup/v_old.dat'//C_NULL_CHAR)
    ren=rename(trim(folder)//'/backup/w.dat'//C_NULL_CHAR,trim(folder)//'/backup/w_old.dat'//C_NULL_CHAR)
#if phiflag == 1
    ren=rename(trim(folder)//'/backup/phi.dat'//C_NULL_CHAR,trim(folder)//'/backup/phi_old.dat'//C_NULL_CHAR)
#if psiflag == 1
    ren=rename(trim(folder)//'/backup/psi_fg.dat'//C_NULL_CHAR,trim(folder)//'/backup/psi_fg_old.dat'//C_NULL_CHAR)
#endif
#endif
#if tempflag == 1
    ren=rename(trim(folder)//'/backup/T.dat'//C_NULL_CHAR,trim(folder)//'/backup/T_old.dat'//C_NULL_CHAR)
#endif
#endif
#if particles == 1
    do j=1,nset
     write(setnum,'(i3.3)') j
     ren=rename(trim(folder)//'/backup/pos_'//setnum//'.dat'//C_NULL_CHAR,trim(folder)//'/backup/pos_'//setnum//'_old.dat'//C_NULL_CHAR)
     ren=rename(trim(folder)//'/backup/vel_'//setnum//'.dat'//C_NULL_CHAR,trim(folder)//'/backup/vel_'//setnum//'_old.dat'//C_NULL_CHAR)
    enddo
#endif
  endif
  in=trim(folder)//'/time_check.dat'
  out=trim(folder)//'/backup/time_check.dat'
  call copy_append(in,out,nt)
#if stat_dump > 0
  ren=rename(trim(folder)//'/backup/stats.dat'//C_NULL_CHAR,trim(folder)//'/backup/stats_old.dat'//C_NULL_CHAR)
  ren=rename(trim(folder)//'/backup/budget.dat'//C_NULL_CHAR,trim(folder)//'/backup/budget_old.dat'//C_NULL_CHAR)
  ren=rename(trim(folder)//'/backup/power_xspectra.dat'//C_NULL_CHAR,trim(folder)//'/backup/power_xspectra_old.dat'//C_NULL_CHAR)
  ren=rename(trim(folder)//'/backup/power_yspectra.dat'//C_NULL_CHAR,trim(folder)//'/backup/power_yspectra_old.dat'//C_NULL_CHAR)
  in=trim(folder)//'/stats.dat'
  out=trim(folder)//'/backup/stats.dat'
  call copy(in,out)
  in=trim(folder)//'/budget.dat'
  out=trim(folder)//'/backup/budget.dat'
  call copy(in,out)
  in=trim(folder)//'/power_xspectra.dat'
  out=trim(folder)//'/backup/power_xspectra.dat'
  call copy(in,out)
  in=trim(folder)//'/power_yspectra.dat'
  out=trim(folder)//'/backup/power_yspectra.dat'
  call copy(in,out)
#endif
! if machine != 4
#else
  if(check.eqv..true.)then
    call rename(trim(folder)//'/backup/time_check.dat',trim(folder)//'/backup/time_check_old.dat')
    call rename(trim(folder)//'/backup/iteration.dat',trim(folder)//'/backup/iteration_old.dat')
#if spectral_flag == 1
    call rename(trim(folder)//'/backup/uc.dat',trim(folder)//'/backup/uc_old.dat')
    call rename(trim(folder)//'/backup/vc.dat',trim(folder)//'/backup/vc_old.dat')
    call rename(trim(folder)//'/backup/wc.dat',trim(folder)//'/backup/wc_old.dat')
#if phiflag == 1
    call rename(trim(folder)//'/backup/phic.dat',trim(folder)//'/backup/phic_old.dat')
#if psiflag == 1
    call rename(trim(folder)//'/backup/psic_fg.dat',trim(folder)//'/backup/psic_fg_old.dat')
#endif
#endif
#if tempflag == 1
    call rename(trim(folder)//'/backup/Tc.dat',trim(folder)//'/backup/Tc_old.dat')
#endif
#elif spectral_flag == 0
    call rename(trim(folder)//'/backup/u.dat',trim(folder)//'/backup/u_old.dat')
    call rename(trim(folder)//'/backup/v.dat',trim(folder)//'/backup/v_old.dat')
    call rename(trim(folder)//'/backup/w.dat',trim(folder)//'/backup/w_old.dat')
#if phiflag == 1
    call rename(trim(folder)//'/backup/phi.dat',trim(folder)//'/backup/phi_old.dat')
#if psiflag == 1
    call rename(trim(folder)//'/backup/psi_fg.dat',trim(folder)//'/backup/psi_fg_old.dat')
#endif
#endif
#if tempflag == 1
    call rename(trim(folder)//'/backup/T.dat',trim(folder)//'/backup/T_old.dat')
#endif
#endif
#if particles == 1
    do j=1,nset
     write(setnum,'(i3.3)') j
     call rename(trim(folder)//'/backup/pos_'//setnum//'.dat',trim(folder)//'/backup/pos_'//setnum//'_old.dat')
     call rename(trim(folder)//'/backup/vel_'//setnum//'.dat',trim(folder)//'/backup/vel_'//setnum//'_old.dat')
    enddo
#endif
  endif
  call system('cp '//trim(folder)//'/time_check.dat '//trim(folder)//'/backup/time_check.dat')
#if stat_dump > 0
  call rename(trim(folder)//'/backup/stats.dat',trim(folder)//'/backup/stats_old.dat')
  call rename(trim(folder)//'/backup/budget.dat',trim(folder)//'/backup/budget_old.dat')
  call rename(trim(folder)//'/backup/power_xspectra.dat',trim(folder)//'/backup/power_xspectra_old.dat')
  call rename(trim(folder)//'/backup/power_yspectra.dat',trim(folder)//'/backup/power_yspectra_old.dat')
  call system('cp '//trim(folder)//'/stats.dat '//trim(folder)//'/backup/stats.dat')
  call system('cp '//trim(folder)//'/budget.dat '//trim(folder)//'/backup/budget.dat')
  call system('cp '//trim(folder)//'/power_xspectra.dat '//trim(folder)//'/backup/power_xspectra.dat')
  call system('cp '//trim(folder)//'/power_yspectra.dat '//trim(folder)//'/backup/power_yspectra.dat')
#endif
#endif
  open(5,file=trim(folder)//'/backup/iteration.dat',status='unknown',form='formatted')
  write(5,'(a,i8)') 'Checkpoint at iteration ',nt
  close(5,status='keep')
endif


if(rank.lt.flow_comm_lim)then
#if spectral_flag == 1
  call write_output_recovery_s(uc,'uc   ')
  call write_output_recovery_s(vc,'vc   ')
  call write_output_recovery_s(wc,'wc   ')
#if phiflag == 1
  call write_output_recovery_s(phic,'phic ')
#if psiflag == 1
  call fine2coarse(psic_fg,psic)
  call write_output_recovery_s_fg(psic_fg,'psic ')
#endif
#endif
#if tempflag == 1
  call write_output_recovery_s(thetac,'Tc   ')
#endif
#elif spectral_flag == 0
  call spectral_to_phys(uc,u,0)
  call spectral_to_phys(vc,v,0)
  call spectral_to_phys(wc,w,0)
  call write_output_recovery(u,'u    ')
  call write_output_recovery(v,'v    ')
  call write_output_recovery(w,'w    ')
#if phiflag == 1
  call spectral_to_phys(phic,phi,0)
  call write_output_recovery(phi,'phi  ')
#if psiflag == 1
  call fine2coarse(psic_fg,psic)
  call spectral_to_phys(psic,psi,0)
  call write_output_recovery_fg(psi_fg,'psi  ')
#endif
#endif
#if tempflag == 1
  call spectral_to_phys(thetac,theta,0)
  call write_output_recovery(theta,'T    ')
#endif
#endif
endif

#if particles == 1
if(rank.ge.leader)then
  call write_output_part_recovery(xp,'pos  ')
  call write_output_part_recovery(up,'vel  ')
endif
#endif

if(rank.eq.0) write(*,*) 'Checkpoint written'

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_recovery(u,namevar)

use mpi
use commondata
use par_size
use mpiIo

double precision :: u(nx,fpz,fpy)

integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=40) :: fname
character(len=5) :: namevar


fname=trim(folder)//'/backup/'//trim(namevar)//'.dat'

offset=0

call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

call mpi_file_write_all(f_handle,u,nx*fpy*fpz,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_recovery_s(uc,namevar)

use mpi
use commondata
use par_size
use mpiIo
use shrink_grid

double precision :: uc(spx,nz,spy,2),ucd(dimc(1),dimc(2),dimc(3),2)

integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=40) :: fname
character(len=5) :: namevar


call shrink_domain(uc,ucd)

fname=trim(folder)//'/backup/'//trim(namevar)//'.dat'

offset=0

! call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
call mpi_file_open(sp_save_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

! call mpi_file_write_all(f_handle,uc,spx*spy*nz*2,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_write_all(f_handle,ucd,dimc(1)*dimc(2)*dimc(3)*2,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_recovery_fg(u,namevar)

use mpi
use commondata
use par_size
use mpiIo
use dual_grid

double precision :: u(npsix,fpzpsi,fpypsi)

integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=40) :: fname
character(len=5) :: namevar


fname=trim(folder)//'/backup/'//trim(namevar)//'_fg.dat'

offset=0

call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype_fg,'native',mpi_info_null,ierr)

call mpi_file_write_all(f_handle,u,npsix*fpypsi*fpzpsi,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_recovery_s_fg(uc,namevar)

use mpi
use commondata
use par_size
use mpiIo
use dual_grid
use shrink_grid

double precision :: uc(spxpsi,npsiz,spypsi,2),ucd(dimc_fg(1),dimc_fg(2),dimc_fg(3),2)

integer :: f_handle ! file handle

integer(mpi_offset_kind) :: offset

character(len=40) :: fname
character(len=5) :: namevar


call shrink_domain_fg(uc,ucd)

fname=trim(folder)//'/backup/'//trim(namevar)//'_fg.dat'

offset=0

! call mpi_file_open(flow_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)
call mpi_file_open(sp_save_comm_fg,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)


!call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'external32',mpi_info_null,ierr)
!call mpi_file_set_view(f_handle,offset,mpi_double_precision,stype_fg,'internal',mpi_info_null,ierr)
call mpi_file_set_view(f_handle,offset,mpi_double_precision,ftype,'native',mpi_info_null,ierr)

! call mpi_file_write_all(f_handle,uc,spxpsi*spypsi*npsiz*2,mpi_double_precision,mpi_status_ignore,ierr)
call mpi_file_write_all(f_handle,ucd,dimc_fg(1)*dimc_fg(2)*dimc_fg(3)*2,mpi_double_precision,mpi_status_ignore,ierr)


call mpi_file_close(f_handle,ierr)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_output_part_recovery(var,namevar)

use mpi
use mpiIO
use commondata
use particle

integer :: j
integer :: f_handle ! file handle
integer(mpi_offset_kind) :: offset
character(len=8) :: time
character(len=40) :: fname
character(len=5) :: namevar
character(len=3) :: setnum

double precision :: var(part_number,3,nset),varloc(part_index(rank_loc+1,2),3)

do j=1,nset
  write(setnum,'(i3.3)') j

  fname=trim(folder)//'/backup/'//trim(namevar)//'_'//setnum//'.dat'

  offset=0
  varloc=var(part_index(rank_loc+1,1)+1:part_index(rank_loc+1,1)+part_index(rank_loc+1,2),:,j)

  call mpi_file_open(part_comm,fname,mpi_mode_create+mpi_mode_rdwr,mpi_info_null,f_handle,ierr)

  !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'external32',mpi_info_null,ierr)
  !call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'internal',mpi_info_null,ierr)
  call mpi_file_set_view(f_handle,offset,mpi_double_precision,part_save,'native',mpi_info_null,ierr)

  call mpi_file_write_all(f_handle,varloc,part_index(rank_loc+1,2)*3,mpi_double_precision,mpi_status_ignore,ierr)

  call mpi_file_close(f_handle,ierr)
enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine copy(in,out)

integer :: outcome

character(len=500) :: buf
character(len=80) :: in,out

logical :: check

inquire(file=trim(in),exist=check)

if(check.eqv..true.)then
  open(12,file=trim(out),status='replace',form='formatted')
  open(11,file=trim(in),status='old',form='formatted')

  do
  read(11,'(a)',iostat=outcome) buf
  if(outcome.eq.0)then
    write(12,'(a)') trim(buf)
  else
    exit
  endif
  enddo

  close(11,status='keep')
  close(12,status='keep')
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine copy_append(in,out,i)

use sim_par

integer :: i,j,outcome

character(len=500) :: buf
character(len=80) :: in,out

logical :: check

inquire(file=trim(in),exist=check)

if(check.eqv..true.)then
  open(12,file=trim(out),status='old',form='formatted',position='append')
  open(11,file=trim(in),status='old',form='formatted')
  do j=1,i-dump_failure+2
    read(11,*)
  enddo

  do j=1,dump_failure
   read(11,'(a)') buf
   write(12,'(a)') trim(buf)
  enddo

  close(11,status='keep')
  close(12,status='keep')
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine copy_append_old(in,out,i)

use sim_par

integer :: i,j,outcome

character(len=500) :: buf
character(len=80) :: in,out

logical :: check

inquire(file=trim(in),exist=check)

if(check.eqv..true.)then
  open(12,file=trim(out),status='old',form='formatted',position='append')
  open(11,file=trim(in),status='old',form='formatted')
  do j=1,i-2*dump_failure+2
    read(11,*)
  enddo

  do j=1,dump_failure
   read(11,'(a)',iostat=outcome) buf
   if(outcome.eq.0) write(12,'(a)') trim(buf)
  enddo

  close(11,status='keep')
  close(12,status='keep')
endif

return
end
