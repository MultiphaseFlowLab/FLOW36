module commondata
 implicit none
 integer :: nx,ny,nz
 integer, parameter :: rank=0
! optional parameter for MPI subroutine with USE MPI_F08, included so that the code is retrocompatible with USE MPI
 integer :: ierr
! for use mpi_f08
!type(mpi_comm) :: cart_comm
 integer :: cart_comm
 double precision :: xl,yl
 character(len=50) :: folder='./output/'
end module commondata



module grid
 use commondata
 double precision, allocatable :: x(:),y(:),z(:)
end module grid



module velocity
 double precision, allocatable, dimension(:,:,:) :: u, v, w
end module velocity



module sim_par
 double precision, parameter :: pi=3.14159265358979
 double precision :: Re, dt, gradpx, gradpy, Co
 double precision :: gamma
! b.c. on w and omega for Helmholtz problem
 double precision, dimension(2) :: p_u,q_u,r_u,p_v,q_v,r_v,p_w,q_w,r_w,p_o,q_o,r_o,zp
 integer :: bc_up, bc_low
 integer :: restart, nt_restart, in_cond, dump_failure
 integer :: nstart, nend, ndump,sdump
end module sim_par



module phase_field
 integer :: phi_flag,in_cond_phi,b_type
 double precision :: rhor,visr,We,Ch,Pe,Fr,grav(3),s_coeff
 double precision, allocatable :: phi(:,:,:), phic(:,:,:,:)
end module phase_field



module surfactant
! mean surfactant concentration, initial condition parameters, ...
 integer :: psi_flag, in_cond_psi
 double precision :: Ex, Pe_psi, P_i, El
 double precision, allocatable :: psi(:,:,:),psic(:,:,:,:)
end module surfactant



module temperature
 integer :: theta_flag,in_cond_theta
 double precision :: Ra,Pr
 double precision, dimension(2) :: p_theta,q_theta,r_theta
 double precision, allocatable :: theta(:,:,:),thetac(:,:,:,:)
end module temperature



module stats
 integer :: flowiter,stat_dump,stat_start
 integer :: plane_comm,col_comm
end module stats
