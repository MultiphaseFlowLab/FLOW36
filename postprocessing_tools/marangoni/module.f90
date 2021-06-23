module commondata
 integer :: nx,ny,nz,nxfg,nyfg,nzfg
 integer :: nstart,nend,delta
 integer :: phi_flag,psi_flag
 integer :: matchedrho,matchedvis
 integer :: expx,expy,expz
 integer :: ierr,ntask,rank
 integer :: spectral

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,ch,we,betas,visr,rhor,dt
 double precision :: xl,yl
 double precision, allocatable, dimension(:) :: x,y,z
end module commondata

module fields
 double precision, dimension(:,:,:), allocatable :: u,v,w,phi,psi
 double precision, dimension(:,:,:,:), allocatable :: uc,vc,wc,phic,psic
end module fields

module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
 type(c_ptr) :: plan_x_fwd_fg,plan_y_fwd_fg,plan_z_fwd_fg
 type(c_ptr) :: plan_x_bwd_fg,plan_y_bwd_fg,plan_z_bwd_fg
 type(c_ptr) :: plan_x_fwd_2D,plan_y_fwd_2D
end module fftw3

module wavenumber
 use commondata
 double precision, allocatable :: kx(:),ky(:),k2(:,:),kxfg(:),kyfg(:),k2fg(:,:)
end module wavenumber
