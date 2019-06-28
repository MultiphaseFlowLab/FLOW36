program FLOW_36

use mpi
use commondata
use mpiIO
use par_size
use grid
use velocity
use sim_par
use phase_field
use stats
use surfactant
use temperature
use particle
use wavenumber
use comm_pattern

#define machine machineflag
#define particles particlecompflag

integer :: i,j,k
integer :: dims(2)
integer :: g_size(3),s_size(3)
integer :: fysize(nycpu),fzsize(nzcpu)
integer :: cysize(nzcpu),cxsize(nycpu)
integer :: l_comm

double precision :: stime,etime,dtime,mtime,gstime,getime,time
double precision :: dx,dy
double precision :: int_1

character(len=5) :: namevar
character(len=50) :: string

logical :: periodic(2), reorder



! start MPI
call mpi_init(ierr)


call split_comm


! check min(nx,ny,nz)>=max(nycpu,nzcpu)
if(min(nx/2,ny).lt.nycpu.or.min(ny,nz).lt.nzcpu) then
  if(rank.eq.0) write(*,*) '#################################'
  if(rank.eq.0) write(*,*) ' STOP, LESS THAN 1 POINT PER CPU'
  if(rank.eq.0) write(*,*) '#################################'
else

  ! read simulation parameters, done in parallel from each task, and print
  ! initial simulation parameters
  call read_input


  ! define array sizes
  call define_sizes

if(rank.lt.flow_comm_lim)then
  ! create cartesian communicator
  dims(1)=nzcpu
  dims(2)=nycpu
  ! for use mpi_f08
  !periodic(1)=1
  !periodic(2)=1
  !reorder=0
  periodic(1)=.true.
  periodic(2)=.true.
  reorder=.false.

  call mpi_cart_create(flow_comm,2,dims,periodic,reorder,cart_comm,ierr)


  ! create derived datatype used in MPI I/O and commit it
  call mpi_cart_sub(cart_comm,[.false.,.true.],l_comm,ierr)
  call mpi_allgather(fpy,1,mpi_integer,fysize,1,mpi_integer,l_comm,ierr)
  call mpi_allgather(spx,1,mpi_integer,cxsize,1,mpi_integer,l_comm,ierr)
  call mpi_comm_free(l_comm,ierr)

  call mpi_cart_sub(cart_comm,[.true.,.false.],l_comm,ierr)
  call mpi_allgather(fpz,1,mpi_integer,fzsize,1,mpi_integer,l_comm,ierr)
  call mpi_allgather(spy,1,mpi_integer,cysize,1,mpi_integer,l_comm,ierr)
  call mpi_comm_free(l_comm,ierr)

  g_size=[nx, nz, ny]
  s_size=[nx, fpz, fpy]
  fstart=[0, 0, 0]
  cstart=[0, 0, 0]
  do i=1,mod(rank,nycpu)
    fstart(3)=fstart(3)+fysize(i)
    cstart(1)=cstart(1)+cxsize(i)
  enddo
  do i=1,floor(real(rank)/real(nycpu))
    fstart(2)=fstart(2)+fzsize(i)
    cstart(3)=cstart(3)+cysize(i)
  enddo
! physical space saving
  call mpi_type_create_subarray(3,g_size,s_size,fstart,mpi_order_fortran, &
   &     mpi_double_precision,ftype,ierr)

  call mpi_type_commit(ftype,ierr)

! spectral space saving (coarse grid)
  call shrink_mapping


! define dual grid sizes
  call define_dual_grid(cxsize,cysize)

! spectral space saving (fine grid)
  call shrink_mapping_fg

! define dual grid communication pattern
  call define_address
endif

  ! create fftw plans
  if(rank.eq.0) then
    write(*,*) 'Starting FFTW plan creation ...'
    write(*,*) 'This may take a while'
  endif
  stime=mpi_wtime()
  call create_plan
  call create_plan_fg
  etime=mpi_wtime()
  dtime=etime-stime
  call mpi_reduce(dtime,mtime,1,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
  mtime=mtime/dble(ntask)
  string='Time elapsed for FFTW plan creation='
  if(rank.eq.0) call write_time(string,mtime)


  !define grid
  if(nx.gt.1.and.ny.gt.1)then
    dx=xl/dble(nx-1)
    dy=yl/dble(ny-1)
  else
    write(*,*) 'error, not enough points'
    stop
  endif

  x(1)=0.0d0
  y(1)=0.0d0
  do i=2,nx
    x(i)=x(i-1)+dx
  enddo
  do j=2,ny
    y(j)=y(j-1)+dy
  enddo

  do k = 1, nz
    z(k)=dcos(((k-1)*pi)/(nz-1))
  enddo


  ! save grid arrays
  if(rank.eq.0) call dump_grid



  ! calculate wavenumbers
  call wave_numbers


! executed only by flow_comm
if(rank.lt.flow_comm_lim)then
  ! initialize flow field according to initial condition provided
  ! allocate velocity matrices in physical and modal space
  call initialize

#define phiflag phicompflag
#define psiflag psicompflag
#if phiflag == 1
  call initialize_phi
#if psiflag == 1
  call initialize_psi
#endif
#endif

#define tempflag tempcompflag
#if tempflag == 1
  call initialize_temp
#endif

#define stat_dump stats_dump_frequency
#if stat_dump > 0
    call initialize_stats
#endif


  ! save initial fields
  if(restart.eq.0)then
    call spectral_to_phys(uc,u,0)
    call spectral_to_phys(vc,v,0)
    call spectral_to_phys(wc,w,0)
    namevar='u'
    call write_output(u,nstart,namevar)
    namevar='v'
    call write_output(v,nstart,namevar)
    namevar='w'
    call write_output(w,nstart,namevar)
#if machine != 2 && machine != 5
    namevar='uc'
    call write_output_spectral(uc,nstart,namevar)
    namevar='vc'
    call write_output_spectral(vc,nstart,namevar)
    namevar='wc'
    call write_output_spectral(wc,nstart,namevar)
#endif
#if phiflag == 1
    call spectral_to_phys(phic,phi,0)
    namevar='phi'
    call write_output(phi,nstart,namevar)
#if machine != 2 && machine != 5
    namevar='phic'
    call write_output_spectral(phic,nstart,namevar)
#endif
#if psiflag == 1
    call spectral_to_phys_fg(psic_fg,psi_fg,0)
    namevar='psi'
    call write_output_fg(psi_fg,nstart,namevar)
#if machine != 2 && machine != 5
    namevar='psic'
    call write_output_spectral_fg(psic_fg,nstart,namevar)
#endif
#endif
#endif
#if tempflag == 1
    call spectral_to_phys(thetac,theta,0)
    namevar='T'
    call write_output(theta,nstart,namevar)
#if machine != 2 && machine != 5
    namevar='Tc'
    call write_output_spectral(thetac,nstart,namevar)
#endif
#endif
  endif

#if phiflag == 1
  call integral_phi(int_1)
#endif

if(restart.eq.0)then
  if(rank.eq.0) call initialize_check(int_1)
endif

endif

if(part_flag.eq.1)then
  ! executed only by part_comm, allocates variables
  if(rank.ge.leader) call allocate_particle
  call create_communication_pattern
  call initialize_particle
endif

  gstime=mpi_wtime()

  ! loop over time
  time=0.0d0
  string='time elapsed per timestep'
  do i=nstart+1,nend

    stime=mpi_wtime()

    time=time+dt


    ! write to screen simulation status
    if(rank.eq.0)then
      write(*,*) '-----------------------------------------------------------------------'
      write(*,'(1x,a,i8,a,i8,2(8x,a,es9.3))') 'step',i,' of ',nend,'t-=',time,'t+=',time*Re
    endif


    call solver(i)


   if(rank.lt.flow_comm_lim)then
#if stat_dump > 0
    if(mod(i,stat_dump).eq.0.and.i.ge.stat_start.and.i.ne.nstart) then
      if(rank.eq.0) write(*,*) 'writing statistics to file'
      call statistics
    endif
#endif

#define fdump physical_dump_frequency
#if fdump > 0
     ! save fields at end of timestep according to ndump value
     if(mod(i,ndump).eq.0.and.i.ne.nstart) then
       if(rank.eq.0) write(*,*) 'saving solution in physical space'
       call spectral_to_phys(uc,u,0)
       namevar='u'
       call write_output(u,i,namevar)
       call spectral_to_phys(vc,v,0)
       namevar='v'
       call write_output(v,i,namevar)
       call spectral_to_phys(wc,w,0)
       namevar='w'
       call write_output(w,i,namevar)
#if phiflag == 1
       call spectral_to_phys(phic,phi,0)
       namevar='phi'
       call write_output(phi,i,namevar)
#if psiflag == 1
       call spectral_to_phys_fg(psic_fg,psi_fg,0)
       namevar='psi'
       call write_output_fg(psi_fg,i,namevar)
#endif
#endif
#if tempflag == 1
       call spectral_to_phys(thetac,theta,0)
       namevar='T'
       call write_output(theta,i,namevar)
#endif
     endif
#endif

#define sdump spectral_dump_frequency
#if sdump > 0
#if machine != 2 && machine != 5
     ! save fields at end of timestep according to ndump value
     if(mod(i,sdump).eq.0.and.i.ne.nstart) then
       if(rank.eq.0) write(*,*) 'saving solution in spectral space'
       namevar='uc'
       call write_output_spectral(uc,i,namevar)
       namevar='vc'
       call write_output_spectral(vc,i,namevar)
       namevar='wc'
       call write_output_spectral(wc,i,namevar)
#if phiflag == 1
       namevar='phic'
       call write_output_spectral(phic,i,namevar)
#if psiflag == 1
       namevar='psic'
       call write_output_spectral_fg(psic_fg,i,namevar)
#endif
#endif
#if tempflag == 1
       namevar='Tc'
       call write_output_spectral(thetac,i,namevar)
#endif
     endif
#endif
#endif

#if phiflag == 1
     call integral_phi(int_1)
#endif
    endif

    if(rank.eq.0) call sim_check(i,int_1)

    ! the MPI all reduce is too costly, so it has been replaced by mean time rank 0
    ! write to screen time elapsed at each time step
    etime=mpi_wtime()
    dtime=etime-stime
!    call mpi_reduce(dtime,mtime,1,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
!    mtime=mtime/dble(ntask)
!    if(rank.eq.0)then
!      call write_time(string,mtime)
!    endif
     if(rank.eq.0) call write_time(string,dtime)


    ! write temporary output files from which the simulation can be restarted if it crashes or stops
    if(mod(i,dump_failure).eq.0)then
      call write_failure(i)
    endif


  enddo

  if(rank.eq.0)then
    write(*,*) '-----------------------------------------------------------------------'
    write(*,*) 'End of simulation'
  endif

  if(rank.lt.flow_comm_lim)then
#if stat_dump > 0
   if(mod(nend,stat_dump).ne.0) then
     call statistics
   endif
   call del_old_stats
#endif

   if(mod(nend,ndump).ne.0.or.ndump.lt.0) then
     if(rank.eq.0) write(*,*) 'Writing final fields in physical space'
     call spectral_to_phys(uc,u,0)
     namevar='u'
     call write_output(u,nend,namevar)
     call spectral_to_phys(vc,v,0)
     namevar='v'
     call write_output(v,nend,namevar)
     call spectral_to_phys(wc,w,0)
     namevar='w'
     call write_output(w,nend,namevar)
#if phiflag == 1
     call spectral_to_phys(phic,phi,0)
     namevar='phi'
     call write_output(phi,nend,namevar)
#if psiflag == 1
     call spectral_to_phys_fg(psic_fg,psi_fg,0)
     namevar='psi'
     call write_output_fg(psi_fg,nend,namevar)
#endif
#endif
#if tempflag == 1
     call spectral_to_phys(thetac,theta,0)
     namevar='T'
     call write_output(theta,nend,namevar)
#endif
   endif

#if machine != 2 && machine != 5
   if(mod(nend,sdump).ne.0.or.sdump.lt.0) then
     if(rank.eq.0) write(*,*) 'Writing final fields in spectral space'
     namevar='uc'
     call write_output_spectral(uc,nend,namevar)
     namevar='vc'
     call write_output_spectral(vc,nend,namevar)
     namevar='wc'
     call write_output_spectral(wc,nend,namevar)
#if phiflag == 1
     namevar='phic'
     call write_output_spectral(phic,nend,namevar)
#if psiflag == 1
     namevar='psic'
     call write_output_spectral_fg(psic_fg,nend,namevar)
#endif
#endif
#if tempflag == 1
     namevar='Tc'
     call write_output_spectral(thetac,nend,namevar)
#endif
   endif
#endif
  endif

  ! output to screen total time elapsed
  getime=mpi_wtime()
  dtime=getime-gstime

  call mpi_reduce(dtime,mtime,1,mpi_double,mpi_sum,0,mpi_comm_world,ierr)
  mtime=mtime/dble(ntask)
  string='total time elapsed for time advancement='
  if(rank.eq.0) call write_time(string,mtime)


  if(rank.lt.flow_comm_lim)then
   ! deallocate all allocatable arrays
   call destroy
#if phiflag == 1
   call destroy_phi
#if psiflag == 1
   call destroy_psi
#endif
#endif
#if tempflag == 1
   call destroy_theta
#endif
  ! destroy MPI I/O derived datatype
   call mpi_type_free(ftype,ierr)
   call mpi_type_free(ftype_fg,ierr)
   if(sp_save_comm.ne.MPI_comm_null) then
     call mpi_type_free(stype,ierr)
     call mpi_comm_free(sp_save_comm,ierr)
   endif
   if(sp_save_comm_fg.ne.MPI_comm_null) then
     call mpi_type_free(stype_fg,ierr)
     call mpi_comm_free(sp_save_comm_fg,ierr)
   endif
   ! destroy cartesian communicator
   call mpi_comm_free(cart_comm,ierr)
  endif
  deallocate(kxpsi)
  deallocate(kypsi)
  deallocate(k2psi)

  ! destroy fftw plans
  call destroy_plan

  ! free all group communicators
  if(rank.lt.flow_comm_lim) call mpi_comm_free(flow_comm,ierr)
#if particles == 1
  if(rank.ge.leader) call mpi_comm_free(part_comm,ierr)
  if(rank.le.flow_comm_lim) call mpi_comm_free(comm_comm,ierr)
  if(rank.ge.leader)then
   call mpi_win_fence(0,window_u,ierr)
   call mpi_win_fence(0,window_v,ierr)
   call mpi_win_fence(0,window_w,ierr)
   call mpi_win_fence(0,window_fx,ierr)
   call mpi_win_fence(0,window_fy,ierr)
   call mpi_win_fence(0,window_fz,ierr)
   call mpi_win_free(window_u,ierr)
   call mpi_win_free(window_v,ierr)
   call mpi_win_free(window_w,ierr)
   call mpi_win_free(window_fx,ierr)
   call mpi_win_free(window_fy,ierr)
   call mpi_win_free(window_fz,ierr)
   deallocate(part_index)
   deallocate(saved_size,address_start)
  endif
  ! call free(xp)
  ! call free(up)
  ! call free(uf)
  ! call free(vf)
  ! call free(wf)
  ! call free(fb_x)
  ! call free(fb_y)
  ! call free(fb_z)
#endif

endif

call mpi_finalize(ierr)


end program FLOW_36
