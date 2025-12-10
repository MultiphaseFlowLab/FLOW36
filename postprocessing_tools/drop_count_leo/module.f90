module commondata
 integer, save :: nx, ny, nz
 integer :: rank, ntask
 integer :: nstart,nend,dump
 integer :: spectral,phiflag,psiflag,b_type
 integer :: counter

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,rhor,visr,we,ch,fr,pe,pe_psi,Ex,P_i,el,grav(3),gradpx,gradpy,dt
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: phi
 double precision, allocatable, dimension(:,:,:,:) :: phic
 double precision, allocatable,dimension(:) :: x,y,z
end module commondata


!module fftw3
! use, intrinsic :: iso_c_binding
! include 'fftw3.f03'
! type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
! type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
!end module fftw3


module velocity
 double precision, allocatable, dimension(:,:,:) :: psi
 double precision, allocatable, dimension(:,:,:,:) :: psic
end module velocity
