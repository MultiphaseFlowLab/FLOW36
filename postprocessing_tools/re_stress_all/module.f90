module commondata 
 integer :: nx,nz,ny
 integer :: st_count
 integer :: count
 double precision :: re
 double precision, dimension(:), allocatable :: greuu,greuv,greuw,grevv,grevw,greww
end module commondata

module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
end module fftw3
