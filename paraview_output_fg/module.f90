module commondata
 integer :: nx, ny, nz
 integer :: nxf,nyf,nzf
 integer :: rank, ntask
 integer :: nstart,nend,ndump,sdump,part_dump,nset
 integer :: x_start,x_end,dnx
 integer :: y_start,y_end,dny
 integer :: z_start,z_end,dnz
 integer :: uflag,vflag,wflag,upflag,vorflag,strflag,topflag,marflag,div2dflag
 integer :: phiflag,psiflag,tempflag,partposflag,partvelflag
 integer :: part_number
 integer :: spectral
 integer :: exp_x,exp_y,exp_z

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt, betas, ch, we
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi,psi,theta,up,vp,wp,omx,omy,omz,strx,stry,strz
 double precision, allocatable, dimension(:,:,:) :: Qtop,marx,mary,marz,div2d
 double precision, allocatable, dimension(:,:,:,:) :: uc,vc,wc,phic,psic,thetac
 double precision, allocatable, dimension(:) :: x,y,z, xfg,yfg,zfg
 double precision, allocatable, dimension(:,:) :: xpar,upar
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
