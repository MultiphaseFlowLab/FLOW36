module commondata
 implicit none
 integer, parameter :: nx=nnx, ny=nny, nz=nnz
 integer, parameter :: nycpu=nnycpu, nzcpu=nnzcpu
 integer :: rank,ntask,ntask_gl,ntask_sh,rank_loc
 integer :: nodes,leader,flow_comm_lim
! optional parameter for MPI subroutine with USE MPI_F08, included so that the code is retrocompatible with USE MPI
 integer :: ierr
! for use mpi_f08
!type(mpi_comm) :: cart_comm
 integer :: cart_comm,flow_comm,part_comm,comm_comm,lacc_rank
 double precision :: xl,yl
 character(len=50) :: folder='./results'
end module commondata



module dual_grid
 use commondata
 integer, parameter :: exp_x=expansionx, exp_y=expansiony, exp_z=expansionz
 integer :: npsix,npsiy,npsiz
 integer :: cg_size(nycpu,nzcpu,5),fg_size(nycpu,nzcpu,5)
 integer :: spxpsi,spypsi,fpypsi,fpzpsi
 integer :: cstartpsi(3),fstartpsi(3)
 integer :: c2fadd(nycpu*nzcpu,2),f2cadd(nycpu*nzcpu,2)
end module dual_grid



module grid
 use commondata
 double precision :: x(nx),y(ny),z(nz)
end module grid



module velocity
 double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
 double precision, allocatable :: uc(:,:,:,:), vc(:,:,:,:), wc(:,:,:,:)
 double precision, allocatable, dimension(:,:,:) :: u_fg, v_fg, w_fg
 double precision, allocatable, dimension(:,:,:,:) :: uc_fg, vc_fg, wc_fg
 double precision, allocatable :: wa2(:,:,:), wa3(:,:,:)
 double precision, allocatable, dimension(:,:,:,:) :: sgradpx,sgradpy
 double precision, allocatable, dimension(:,:,:) :: forx,fory,forz
end module velocity



module wavenumber
 use commondata
 double precision :: kx(nx/2+1),ky(ny),k2(nx/2+1,ny)
 double precision, allocatable :: kxpsi(:),kypsi(:),k2psi(:,:)
end module wavenumber



module sim_par
 double precision, parameter :: pi=3.14159265358979
 double precision :: Re, dt, gradpx, gradpy, Co, repow
 double precision :: gamma
! b.c. on w and omega for Helmholtz problem
 double precision, dimension(2) :: p_u,q_u,r_u,p_v,q_v,r_v,p_w,q_w,r_w,p_o,q_o,r_o,zp
 integer :: bc_up, bc_low
 integer :: cpi_flag
 integer :: restart, nt_restart, in_cond, dump_failure
 integer :: nstart, nend, ndump,sdump
end module sim_par



module phase_field
 integer :: phi_flag,phicor_flag,in_cond_phi,b_type,phiflux_flag,body_flag,sgradp_flag,ele_flag, non_newtonian
 double precision :: rhor,visr,We,Ch,Pe,Fr,grav(3),s_coeff,lamphi,body_c,body_d(3),sgradp_d(3),stuart
 double precision :: muinfmuzero, exp_non_new
 double precision, allocatable :: phi(:,:,:), phic(:,:,:,:)
 double precision, allocatable :: phi_fg(:,:,:), phic_fg(:,:,:,:)
end module phase_field



module surfactant
! mean surfactant concentration, initial condition parameters, ...
 integer :: psi_flag, in_cond_psi
 double precision :: Ex, Pe_psi, P_i, El
 double precision, allocatable :: psi(:,:,:),psic(:,:,:,:)
 double precision, allocatable :: psi_fg(:,:,:),psic_fg(:,:,:,:)
end module surfactant



module temperature
 integer :: theta_flag,in_cond_theta
 double precision :: Ra,Pr
 double precision, dimension(2) :: p_theta,q_theta,r_theta
 double precision, allocatable :: theta(:,:,:),thetac(:,:,:,:)
end module temperature



module particle
 integer :: part_flag,part_number,in_cond_part_pos,in_cond_part_vel,part_dump,nset
 integer :: subiterations
 integer, allocatable, dimension(:,:) :: part_index
 double precision :: dt_part
 double precision, allocatable, dimension(:) :: stokes,dens_part,d_par
 double precision, pointer, dimension(:,:,:) :: xp,up
 double precision, pointer, dimension(:,:,:) :: uf,vf,wf,Tf,fb_x,fb_y,fb_z
 ! mpi shared memory synchronization windows
 integer :: window_u,window_v,window_w,window_T,window_fx,window_fy,window_fz,window_xp,window_up
end module particle



module velocity_old
 ! store velocities at t=n-1, needed for non-linear part of time derivatives for non-matched densities case
 double precision, allocatable :: ucp(:,:,:,:), vcp(:,:,:,:), wcp(:,:,:,:)
end module velocity_old



module mpiIO
 integer :: ftype,stype
 integer :: ftype_fg,stype_fg
 integer :: sp_save_comm,sp_save_comm_fg
 integer :: part_save,part_save_scalar
end module mpiIo



module par_size
! size in physical space: nx,npz,npy
 integer :: fpy,fpz
! size in spectral space:
 integer :: spx,spy
 integer :: cstart(3),fstart(3)
end module par_size



module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
 type(c_ptr) :: plan_x_fwd_fg,plan_y_fwd_fg,plan_z_fwd_fg
 type(c_ptr) :: plan_x_bwd_fg,plan_y_bwd_fg,plan_z_bwd_fg
end module fftw3




module cufftplans
 use, intrinsic :: iso_c_binding
 integer :: cudaplan_x_fwd, cudaplan_y_fwd, cudaplan_z_fwd
 integer :: cudaplan_x_bwd, cudaplan_y_bwd, cudaplan_z_bwd
 integer :: cudaplan_x_fwd_fg, cudaplan_y_fwd_fg, cudaplan_z_fwd_fg
 integer :: cudaplan_x_bwd_fg, cudaplan_y_bwd_fg, cudaplan_z_bwd_fg
 integer :: cudaplan_z_fwd_1d
 integer :: cudaplan_z_bwd_1d
 integer :: gerr
end module cufftplans


module sterms
 double precision, allocatable, dimension(:,:,:,:) :: s1_o,s2_o,s3_o, sphi_o, spsi_o, stheta_o
end module sterms



module stats
 integer :: flowiter,stat_dump,stat_start
 integer :: plane_comm,col_comm
end module stats



module shrink_grid
 integer :: sx1,sx2,sy1,sy2,sy3,sy4
 integer :: sx1_fg,sx2_fg,sy1_fg,sy2_fg,sy3_fg,sy4_fg
 integer, dimension(3) :: dimc,dimc_fg
 integer :: up,up_fg
end module shrink_grid



module comm_pattern
 integer, dimension(4) :: chunk_size
 integer, allocatable, dimension(:,:) ::saved_size,address_start
end module comm_pattern



!only for profiling using NVTX (Nvidia GPU)
!! only for profiling
!module nvtx

!use iso_c_binding
!implicit none

!integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', Z'00ff0000', Z'00ffffff']
!character,private,target :: tempName(256)

!type, bind(C):: nvtxEventAttributes
!   integer(C_INT16_T):: version=1
!   integer(C_INT16_T):: size=48 !
!   integer(C_INT):: category=0
!   integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
!   integer(C_INT):: color
!   integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
!   integer(C_INT):: reserved0
!   integer(C_INT64_T):: payload   ! union uint,int,double
!   integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1
!   type(C_PTR):: message  ! ascii char
! end type
!
! interface nvtxRangePush
!   ! push range with custom label and standard color
!   subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
!   use iso_c_binding
!   character(kind=C_CHAR) :: name(256)
!   end subroutine
!
!   ! push range with custom label and custom color
!   subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
!   use iso_c_binding
!   import:: nvtxEventAttributes
!   type(nvtxEventAttributes):: event
!   end subroutine
! end interface
!
! interface nvtxRangePop
!   subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
!   end subroutine
! end interface
!
! contains
!
! subroutine nvtxStartRange(name,id)
!   character(kind=c_char,len=*) :: name
!   integer, optional:: id
!   type(nvtxEventAttributes):: event
!   character(kind=c_char,len=256) :: trimmed_name
!   integer:: i
!
!   trimmed_name=trim(name)//c_null_char
!
!   ! move scalar trimmed_name into character array tempName
!   do i=1,LEN(trim(name)) + 1
!      tempName(i) = trimmed_name(i:i)
!   enddo
!
!   if ( .not. present(id)) then
!      call nvtxRangePush(tempName)
!    else
!      event%color=col(mod(id,7)+1)
!      event%message=c_loc(tempName)
!      call nvtxRangePushEx(event)
!    end if
!  end subroutine
!
!  subroutine nvtxEndRange
!    call nvtxRangePop
!  end subroutine
!
!  ! push range with custom label and custom color
!  subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
!  use iso_c_binding
!  import:: nvtxEventAttributes
!  type(nvtxEventAttributes):: event
!  end subroutine
! end interface
!
! interface nvtxRangePop
!  subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
!  end subroutine
! end interface
!
! contains
!
! subroutine nvtxStartRange(name,id)
!  character(kind=c_char,len=*) :: name
!  integer, optional:: id
!  type(nvtxEventAttributes):: event
!  character(kind=c_char,len=256) :: trimmed_name
!  integer:: i
!
!  trimmed_name=trim(name)//c_null_char
!
!  ! move scalar trimmed_name into character array tempName
!  do i=1,LEN(trim(name)) + 1
!     tempName(i) = trimmed_name(i:i)
!  enddo
!
!
!  if ( .not. present(id)) then
!    call nvtxRangePush(tempName)
!  else
!    event%color=col(mod(id,7)+1)
!    event%message=c_loc(tempName)
!    call nvtxRangePushEx(event)
!  end if
! end subroutine
!
! subroutine nvtxEndRange
!  call nvtxRangePop
! end subroutine
!
! end module nvtx
