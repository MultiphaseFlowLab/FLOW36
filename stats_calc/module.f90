module commondata
 integer :: nx, ny, nz
 integer :: rank, ntask
 integer :: nstart,nend,ndump,sdump
 integer :: phiflag
 integer :: spectral
 integer :: counter

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: u
 double precision, allocatable, dimension(:,:,:,:) :: uc
 double precision, allocatable, dimension(:) :: mean_u, rms_u, skw_u, flt_u
 double precision, allocatable,dimension(:) :: x,y,z
end module commondata


module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
end module fftw3


