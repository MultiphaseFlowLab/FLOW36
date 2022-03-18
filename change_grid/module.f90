module commondata
 integer :: nx, ny, nz
 integer :: phiflag,psiflag

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt
end module commondata


module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
end module fftw3


module input
 integer :: inx,iny,inz
end module input


module output
 integer :: onx,ony,onz
end module output

