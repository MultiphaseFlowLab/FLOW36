subroutine hist_term_temp(htheta)

use commondata
use sim_par
use par_size
use temperature
use wavenumber

double precision :: htheta(spx,nz,spy,2)
double precision :: theta_z(spx,nz,spy,2),theta_2z(spx,nz,spy,2)

integer :: i,j

call dz(thetac,theta_z)
call dz(theta_z,theta_2z)

!$acc kernels
do j=1,spy
  do i=1,spx
    ! Crank-Nicolson
    ! htheta(i,:,j,1)=-Re*Pr/dt*htheta(i,:,j,1)-theta_2z(i,:,j,1)-(Re*Pr/dt-k2(i+cstart(1),j+cstart(3)))*thetac(i,:,j,1)
    ! htheta(i,:,j,2)=-Re*Pr/dt*htheta(i,:,j,2)-theta_2z(i,:,j,2)-(Re*Pr/dt-k2(i+cstart(1),j+cstart(3)))*thetac(i,:,j,2)
    ! Implicit Euler
    htheta(i,:,j,1)=-Re*Pr/dt*htheta(i,:,j,1)-Re*Pr/dt*thetac(i,:,j,1)
    htheta(i,:,j,2)=-Re*Pr/dt*htheta(i,:,j,2)-Re*Pr/dt*thetac(i,:,j,2)
  enddo
enddo
!$acc end kernels

return
end
