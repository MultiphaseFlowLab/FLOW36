module commondata
 integer :: nx, ny, nz
 integer :: nxf,nyf,nzf
 integer :: rank, ntask
 integer :: nstart,nend,ndump,sdump
 integer :: x_start,x_end,dnx
 integer :: y_start,y_end,dny
 integer :: z_start,z_end,dnz
 integer :: phiflag, psiflag,tempflag
 integer :: spectral
 integer :: exp_x,exp_y,exp_z

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi,psi,theta
 double precision, allocatable, dimension(:,:,:,:) :: uc,vc,wc,phic,psic,thetac
 double precision, allocatable,dimension(:) :: x,y,z, xfg,yfg,zfg
end module commondata


module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
 type(c_ptr) :: plan_x_fwd_fg,plan_y_fwd_fg,plan_z_fwd_fg
 type(c_ptr) :: plan_x_bwd_fg,plan_y_bwd_fg,plan_z_bwd_fg
end module fftw3
