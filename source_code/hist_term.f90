subroutine hist_term(h1,h2,h3,h)

use commondata
use par_size
use velocity
use wavenumber
use sim_par

double precision, dimension(spx,nz,spy,2) :: h1, h2, h3, h
double precision, dimension(spx,nz,spy,2) :: du,ddu

integer :: i,j


! h1
call dz(uc,du)
call dz(du,ddu)
!$acc parallel loop collapse(2)
do j=1,spy
 do i=1,spx
  h1(i,:,j,1)=h1(i,:,j,1)+gamma*ddu(i,:,j,1)+(1.0d0-k2(i+cstart(1),j+cstart(3))*gamma)*uc(i,:,j,1)
  h1(i,:,j,2)=h1(i,:,j,2)+gamma*ddu(i,:,j,2)+(1.0d0-k2(i+cstart(1),j+cstart(3))*gamma)*uc(i,:,j,2)
 enddo 
enddo 

! h2
call dz(vc,du)
call dz(du,ddu)
!$acc parallel loop collapse(2)
do j=1,spy
 do i=1,spx
  h2(i,:,j,1)=h2(i,:,j,1)+gamma*ddu(i,:,j,1)+(1.0d0-k2(i+cstart(1),j+cstart(3))*gamma)*vc(i,:,j,1)
  h2(i,:,j,2)=h2(i,:,j,2)+gamma*ddu(i,:,j,2)+(1.0d0-k2(i+cstart(1),j+cstart(3))*gamma)*vc(i,:,j,2)
 enddo
enddo

! h3
call dz(wc,du)
call dz(du,ddu)
!$acc parallel loop collapse(2)
do j=1,spy
 do i=1,spx
  h3(i,:,j,1)=h3(i,:,j,1)+gamma*ddu(i,:,j,1)+(1.0d0-k2(i+cstart(1),j+cstart(3))*gamma)*wc(i,:,j,1)
  h3(i,:,j,2)=h3(i,:,j,2)+gamma*ddu(i,:,j,2)+(1.0d0-k2(i+cstart(1),j+cstart(3))*gamma)*wc(i,:,j,2)
 enddo
enddo

! h
!$acc kernels
do j=1,spy
 do i=1,spx
  ddu(i,:,j,1)=-kx(i+cstart(1))*h1(i,:,j,2)-ky(j+cstart(3))*h2(i,:,j,2)
  ddu(i,:,j,2)=kx(i+cstart(1))*h1(i,:,j,1)+ky(j+cstart(3))*h2(i,:,j,1)
 enddo
enddo
!$acc end kernels
call dz(ddu,h)
!$acc kernels
do j=1,spy
 do i=1,spx
  h(i,:,j,1)=h(i,:,j,1)+k2(i+cstart(1),j+cstart(3))*h3(i,:,j,1)
  h(i,:,j,2)=h(i,:,j,2)+k2(i+cstart(1),j+cstart(3))*h3(i,:,j,2)
 enddo
enddo
!$acc end kernels

return
end
