subroutine fine2coarse(varf,varc)

use commondata

double precision :: varc(nx/2+1,nz,ny,2),varf((expx*nx)/2+1,expz*(nz-1)+1,expy*ny,2)

integer :: nfy

nfy=expy*ny

varc=0.0d0

varc(1:nx/2+1,1:nz,1:ny/2+1,1:2)=varf(1:nx/2+1,1:nz,1:ny/2+1,1:2)
varc(1:nx/2+1,1:nz,ny-(ny/2-1)+1:ny,1:2)=varf(1:nx/2+1,1:nz,nfy-(ny/2-1)+1:nfy,1:2)

! renormalize
varc=varc/dble(expx*expy)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coarse2fine(varc,varf)

use commondata

double precision :: varc(nx/2+1,nz,ny,2),varf((expx*nx)/2+1,expz*(nz-1)+1,expy*ny,2)

integer :: nfy

nfy=expy*ny

varf=0.0d0

varf(1:nx/2+1,1:nz,1:ny/2+1,1:2)=varc(1:nx/2+1,1:nz,1:ny/2+1,1:2)
varf(1:nx/2+1,1:nz,nfy-(ny/2-1)+1:nfy,1:2)=varc(1:nx/2+1,1:nz,ny-(ny/2-1)+1:ny,1:2)

! renormalize
varf=varf*dble(expx*expy)

return
end
