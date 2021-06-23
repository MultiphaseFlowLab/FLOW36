module commondata
 integer :: nx, ny, nz
 integer :: nxf,nyf,nzf
 integer :: rank, ntask
 integer :: nstart,nend,ndump,sdump
 integer :: psiflag
 integer :: spectral
 integer :: exp_x,exp_y,exp_z
 integer :: samples_in,samples_out

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt, betas, ch, we
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi,psi
 double precision, allocatable, dimension(:,:,:,:) :: uc,vc,wc,phic,psic
 double precision, allocatable,dimension(:) :: x,y,z
end module commondata


module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
 type(c_ptr) :: plan_x_fwd_fg,plan_y_fwd_fg,plan_z_fwd_fg
 type(c_ptr) :: plan_x_bwd_fg,plan_y_bwd_fg,plan_z_bwd_fg
end module fftw3



module wavenumber
 use commondata
 double precision, allocatable :: kx(:),ky(:)
end module wavenumber



module jpdf
 integer :: nbins
 double precision, allocatable, dimension(:,:,:) :: QvRv,QsRs,QsQo
 double precision, allocatable, dimension(:) :: Qva,Rva,Qsa,Rsa,Qoa
end module jpdf
