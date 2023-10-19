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
 double precision, allocatable, dimension(:,:,:) :: phi,theta,theta_fg
 double precision, allocatable, dimension(:) :: x,y,z, xfg,yfg,zfg, mean
 ! bin for the mean 
 integer :: nbin=200
end module commondata



