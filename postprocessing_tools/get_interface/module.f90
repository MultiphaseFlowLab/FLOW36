module commondata
 integer :: nx, ny, nz
 integer :: rank, ntask
 integer :: nstart,nend,ndump,sdump
 integer :: phiflag
 integer :: spectral
 integer :: counter

 double precision, parameter :: pi=3.14159265358979
 double precision :: re,dt
 double precision :: xl,yl
 double precision, allocatable, dimension(:,:,:) :: phi,w
 double precision, allocatable,dimension(:) :: x,y,z
end module commondata


