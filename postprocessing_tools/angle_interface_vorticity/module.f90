module mpiIO
 integer :: ftype,stype
end module mpiIo


module fftw3
 use, intrinsic :: iso_c_binding
 include 'fftw3.f03'
 type(c_ptr) :: plan_x_fwd,plan_y_fwd,plan_z_fwd
 type(c_ptr) :: plan_x_bwd,plan_y_bwd,plan_z_bwd
end module fftw3


module commondata
 double precision, allocatable :: x(:),y(:),z(:) !,zcm(:)
 double precision :: ucm,vcm,wcm,zcm
 double precision, parameter :: pi=3.14159265358979, threshold=2.0e-1
 integer :: rank,ntask,ierr,cart_comm
 integer :: nycpu,nzcpu,fhandle,total
 integer :: nx,ny,nz,fpy,fpz,spx,spy
 integer, dimension(3) :: fstart,cstart
 integer :: out_spectral,generate_paraview,oldstat
 double precision :: xl,yl
 character(len=80) :: folder
end module commondata


module sim_parameter
 double precision :: Re,Ch,Pe,We,Fr,rhor,visr,grav(3)
 double precision :: Pe_psi,Ex,P_i,El
 integer :: nstart,nend,dump,begin
 integer :: spectral,phiflag,psiflag,b_type
end module sim_parameter


module vars
 double precision, allocatable, dimension(:) :: um,vm,wm
 double precision, allocatable, dimension(:,:,:) :: u,v,w,phi,tke
 double precision, allocatable, dimension(:,:,:,:) :: uc,vc,wc,phic,tkec
end module vars


module wavenumbers
 double precision, allocatable :: kx(:),ky(:),k2(:,:)
end module wavenumbers


module paraview_utils
 integer :: x_start,x_end,dnx
 integer :: y_start,y_end,dny
 integer :: z_start,z_end,dnz
 double precision, allocatable, dimension(:,:,:) :: uvel,vvel,wvel,phi_ph,tkeen
end module paraview_utils



module pdf_calc
 integer, parameter :: n_samples=500
 double precision :: pdf_p(n_samples,2),pdf_t(n_samples,2)
end module pdf_calc
