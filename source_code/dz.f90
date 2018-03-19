subroutine dz(u,duz)
! be careful, do not use same input and output in the call,
! e.g. call dz(u,u) <---- NO!

use commondata
use par_size

double precision :: u(spx,nz,spy,2),duz(spx,nz,spy,2)

integer :: k,m

do m=1,2
  duz(:,nz,:,m)=0.0d0
  duz(:,nz-1,:,m)= 2.0d0*(nz-1)*u(:,nz,:,m)
  do k=nz-2,1,-1
    duz(:,k,:,m)=duz(:,k+2,:,m)+2.0d0*dble(k)*u(:,k+1,:,m)
  enddo
enddo


return
end



subroutine dz_red(u,duz)

use commondata
use par_size

double precision :: u(spx,nz,spy),duz(spx,nz,spy)

integer :: k


duz(:,nz,:)=0.0d0
duz(:,nz-1,:)= 2.0d0*(nz-1)*u(:,nz,:)
do k=nz-2,1,-1
  duz(:,k,:)=duz(:,k+2,:)+2.0d0*dble(k)*u(:,k+1,:)
enddo


return
end
