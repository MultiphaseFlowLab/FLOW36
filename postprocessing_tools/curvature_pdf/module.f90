module commondata
 integer, save :: nx, ny, nz
 integer :: rank, ntask
 integer :: nstart,nend,ndump
 integer :: spectral

 double precision, parameter :: pi=3.14159265358979
 double precision :: re
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: phi,kv
 double precision, allocatable, dimension(:,:,:,:) :: phic
 double precision, allocatable,dimension(:) :: x,y,z
end module commondata


module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
end module fftw3


module wavenumbers
 double precision, allocatable :: kx(:),ky(:),k2(:,:)
end module wavenumbers



module pdf_calc
 double precision, allocatable :: axis(:)
 double precision :: threshold,maxk,mink
 integer, allocatable :: pdf(:)
 integer :: nset
end module pdf_calc
