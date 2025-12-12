module commondata
 integer, save :: nx, ny, nz
 integer :: rank, ntask
 integer :: nstart,nend,dump
 integer :: spectral,phiflag,psiflag,b_type
 integer :: counter

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,rhor,visr,we,ch,fr,pe,pe_psi,Ex,P_i,el,grav(3),gradpx,gradpy,dt
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: phi,phix,phiy,phiz!,div
 double precision, allocatable, dimension(:,:,:,:) :: phic
 double precision, allocatable,dimension(:) :: x,y,z
end module commondata



module velocity
 double precision, allocatable, dimension(:,:,:) :: u,v,w,press,psi
 double precision, allocatable, dimension(:,:,:,:) :: uc,vc,wc,pressc,psic
end module velocity
