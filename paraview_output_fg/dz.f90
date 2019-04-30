subroutine dz(varin,duz)
! be careful, do not use same input and output in the call,
! e.g. call dz(u,u) <---- NO!

use commondata

double precision :: varin(nxf/2+1,nzf,nyf,2),duz(nxf/2+1,nzf,nyf,2)

integer :: k,m

do m=1,2
  duz(:,nzf,:,m)=0.0d0
  duz(:,nzf-1,:,m)= 2.0d0*(nzf-1)*varin(:,nzf,:,m)
  do k=nzf-2,1,-1
    duz(:,k,:,m)=duz(:,k+2,:,m)+2.0d0*dble(k)*varin(:,k+1,:,m)
  enddo
enddo

return
end
