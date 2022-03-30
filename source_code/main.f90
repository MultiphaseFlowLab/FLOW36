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
#define openaccflag openacccompflag
#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define stat_dump stats_dump_frequency
#define particles particlecompflag
#define fdump physical_dump_frequency
#define sdump spectral_dump_frequency
#define cpiflag cpicompflag

#if openaccflag == 1
use openacc
#endif


integer :: i,j,k
integer :: dims(2)
integer :: g_size(3),s_size(3)
integer :: fysize(nycpu),fzsize(nzcpu)
integer :: cysize(nzcpu),cxsize(nycpu)
integer :: l_comm
integer :: numdevices, devicenum

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

! assign GPU to rank
#if openaccflag == 1
  numdevices=acc_get_num_devices(acc_device_nvidia)
  devicenum=mod(rank,numdevices)
  call acc_set_device_num(devicenum,acc_device_nvidia)
! debug only (to be removed)
!  write(*,*) "MPI Rank, GPU number", rank, devicenum
#endif
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


! executed only by flow_comm
if(rank.lt.flow_comm_lim)then
  ! calculate wavenumbers
  call wave_numbers
  ! initialize flow field according to initial condition provided
  ! allocate velocity matrices in physical and modal space
  call initialize

#if phiflag == 1
  call initialize_phi
#if psiflag == 1
  call initialize_psi
#endif
#endif

#if tempflag == 1
  call initialize_temp
#endif

#if stat_dump > 0
    call initialize_stats
#endif

  ! save initial fields
  if(restart.eq.0)then
   call save_flow_comm(nstart)
  endif

#if phiflag == 1
  call integral_phi(int_1)
#endif

 if(restart.eq.0)then
  if(rank.eq.0) call initialize_check(int_1)
 endif
! end of flow_comm only
endif

#if particles == 1
  ! executed only by part_comm, allocates variables
  if(rank.ge.leader) call allocate_particle
  call create_communication_pattern
  call initialize_particle
  ! save particle data (part_comm)
  if(rank.ge.leader)then
   namevar='pos'
   call write_output_part(xp,nstart,namevar)
   namevar='vel'
   call write_output_part(up,nstart,namevar)
   call write_output_partf(nstart)
  endif
#if tempflag == 1
    ! get temperature at particle position (done here, after defining all auxiliary arrays)
    call get_temperature
    if(rank.ge.leader) call write_output_partT(nstart)
#endif
#endif

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

    ! save variables and stats (flow_comm)
    if(rank.lt.flow_comm_lim)then
#if stat_dump > 0
     if(mod(i,stat_dump).eq.0.and.i.ge.stat_start.and.i.ne.nstart) then
      if(rank.eq.0) write(*,*) 'writing statistics to file'
      call statistics
     endif
#endif
     call save_flow_comm(i)
    endif

#if particles == 1
    if(((mod(i,ndump).eq.0.and.ndump.gt.0).or. &
     &      (mod(i,sdump).eq.0.and.sdump.gt.0).or. &
     &      (mod(i,part_dump).eq.0)).and.(i.ne.nstart))then
#if tempflag == 1
     ! get temperature at particle position
     ! done here if particle are temperature tracers (do not fetch temperature at each time step)
     call get_temperature
#endif
     ! save particle data (part_comm)
     if(rank.ge.leader)then
      if(rank.eq.leader) write(*,*) 'saving particle solution'
      namevar='pos'
      call write_output_part(xp,i,namevar)
      namevar='vel'
      call write_output_part(up,i,namevar)
      ! write fluid velocity at particle position
      call write_output_partf(i)
      ! write fluid temperature at particle position
#if tempflag == 1
      call write_output_partT(i)
#endif
     endif
    endif
#endif

#if phiflag == 1
    if(rank.lt.flow_comm_lim) call integral_phi(int_1)
#endif

    if(rank.eq.0) call sim_check(i,int_1)

! only for CPI
#if cpiflag == 1
    if(rank.lt.flow_comm_lim) call mpi_bcast(gradpx,1,mpi_double_precision,0,flow_comm,ierr)
#endif


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

  ! write final fields if not already written
  if(rank.lt.flow_comm_lim)then
#if stat_dump > 0
   if(mod(nend,stat_dump).ne.0) then
     call statistics
   endif
   call del_old_stats
#endif
   call save_flow_comm_final(nend)
  endif

  ! save final particle data if not already written
#if particles == 1
#if tempflag == 1
    ! get temperature at particle position
    call get_temperature
#endif
  if((rank.ge.leader).and. &
 &     ((mod(nend,ndump).ne.0.and.ndump.gt.0).or. &
 &      (mod(nend,sdump).ne.0.and.sdump.gt.0).or. &
 &      (mod(nend,part_dump).ne.0)))then
     namevar='pos'
     call write_output_part(xp,nend,namevar)
     namevar='vel'
     call write_output_part(up,nend,namevar)
     call write_output_partf(nend)
#if tempflag == 1
     call write_output_partT(nend)
#endif
    endif
#endif

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
   deallocate(kxpsi)
   deallocate(kypsi)
   deallocate(k2psi)
  endif


  ! destroy fftw plans
  call destroy_plan

  ! free all group communicators
  if(rank.lt.flow_comm_lim) call mpi_comm_free(flow_comm,ierr)
  deallocate(stokes,dens_part,d_par)
#if particles == 1
  if(rank.le.flow_comm_lim) call mpi_comm_free(comm_comm,ierr)
  if(rank.ge.leader)then
   call mpi_comm_free(part_comm,ierr)
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
   call mpi_type_free(part_save,ierr)
   call mpi_type_free(part_save_scalar,ierr)
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
