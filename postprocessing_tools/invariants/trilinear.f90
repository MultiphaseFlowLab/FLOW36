subroutine trilinear(f,i,j,k,xp,yp,zp,fp)

use commondata

double precision :: f(nx,nz,ny),fp,xp,yp,zp
integer :: i,j,k,ip,jp,kp

double precision :: f00,f01,f10,f11,f0,f1
double precision :: xd,yd,zd


! consider periodicity
xd=(xp-xx(i))/(xx(i+1)-xx(i))
yd=(yp-yy(j))/(yy(j+1)-yy(j))
zd=(zp-zz(k))/(zz(k+1)-zz(k))

ip=i+1
if(ip.gt.nx) ip=ip-nx
jp=j+1
if(jp.gt.ny) jp=jp-ny
kp=k+1
if(kp.gt.nz) return


f00=f(i,k,j)*(1.0d0-xd)+f(ip,k,j)*xd
f01=f(i,k,jp)*(1.0d0-xd)+f(ip,k,jp)*xd
f10=f(i,kp,j)*(1.0d0-xd)+f(ip,kp,j)*xd
f11=f(i,kp,jp)*(1.0d0-xd)+f(ip,kp,jp)*xd

f0=f00*(1.0d0-yd)+f01*yd
f1=f10*(1.0d0-yd)+f11*yd

fp=f0*(1.0d0-zd)+f1*zd

return
end
