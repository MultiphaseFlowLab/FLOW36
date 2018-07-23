subroutine fine2coarse(varf,varc)

use commondata

double precision :: varc(nx/2+1,nz,ny,2),varf((exp_x*nx)/2+1,exp_z*(nz-1)+1,exp_y*ny,2)

integer :: nfy

nfy=exp_y*ny

varc=0.0d0

varc(1:nx/2+1,1:nz,1:ny/2+1,1:2)=varf(1:nx/2+1,1:nz,1:ny/2+1,1:2)
varc(1:nx/2+1,1:nz,ny-(ny/2-1)+1:ny,1:2)=varf(1:nx/2+1,1:nz,nfy-(ny/2-1)+1:nfy,1:2)

! renormalize
varc=varc/dble(exp_x*exp_y)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coarse2fine(varc,varf)

use commondata

double precision :: varc(nx/2+1,nz,ny,2),varf((exp_x*nx)/2+1,exp_z*(nz-1)+1,exp_y*ny,2)

integer :: nfy

nfy=exp_y*ny

varf=0.0d0

varf(1:nx/2+1,1:nz,1:ny/2+1,1:2)=varc(1:nx/2+1,1:nz,1:ny/2+1,1:2)
varf(1:nx/2+1,1:nz,nfy-(ny/2-1)+1:nfy,1:2)=varc(1:nx/2+1,1:nz,ny-(ny/2-1)+1:ny,1:2)

! renormalize
varf=varf*dble(exp_x*exp_y)

return
end
