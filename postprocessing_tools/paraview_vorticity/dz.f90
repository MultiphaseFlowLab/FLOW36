subroutine dz(var,duz)
! be careful, do not use same input and output in the call,
! e.g. call dz(u,u) <---- NO!

use commondata


double precision, dimension(nx/2+1,nz,ny,2) :: var,duz

integer :: k,m

do m=1,2
  duz(:,nz,:,m)=0.0d0
  duz(:,nz-1,:,m)= 2.0d0*(nz-1)*var(:,nz,:,m)
  do k=nz-2,1,-1
    duz(:,k,:,m)=duz(:,k+2,:,m)+2.0d0*dble(k)*var(:,k+1,:,m)
  enddo
enddo


return
end
