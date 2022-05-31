subroutine euler(s1,s2,s3,h1,h2,h3)

use commondata
use sim_par
use par_size

double precision :: s1(spx,nz,spy,2), s2(spx,nz,spy,2), s3(spx,nz,spy,2)
double precision :: h1(spx,nz,spy,2), h2(spx,nz,spy,2), h3(spx,nz,spy,2)

!$acc kernels
h1=dt*s1
h2=dt*s2
h3=dt*s3
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth(s1,s2,s3,h1,h2,h3)

use commondata
use sim_par
use par_size
use sterms

double precision :: s1(spx,nz,spy,2), s2(spx,nz,spy,2), s3(spx,nz,spy,2)
double precision :: h1(spx,nz,spy,2), h2(spx,nz,spy,2), h3(spx,nz,spy,2)

!$acc kernels 
h1=0.5d0*dt*(3.0d0*s1-s1_o)
h2=0.5d0*dt*(3.0d0*s2-s2_o)
h3=0.5d0*dt*(3.0d0*s3-s3_o)
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_phi(sphi,hphi)

use commondata
use sim_par
use par_size

double precision :: sphi(spx,nz,spy,2), hphi(spx,nz,spy,2)

hphi=dt*sphi

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth_phi(sphi,hphi)

use commondata
use sim_par
use par_size
use sterms

double precision :: sphi(spx,nz,spy,2), hphi(spx,nz,spy,2)

hphi=0.5d0*dt*(3.0d0*sphi-sphi_o)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_psi(spsi,hpsi)

use commondata
use sim_par
use par_size
use dual_grid

double precision, dimension(spxpsi,npsiz,spypsi,2) :: spsi, hpsi

hpsi=dt*spsi

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth_psi(spsi,hpsi)

use commondata
use sim_par
use par_size
use sterms
use dual_grid

double precision, dimension(spxpsi,npsiz,spypsi,2) :: spsi, hpsi

hpsi=0.5d0*dt*(3.0d0*spsi-spsi_o)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_theta(stheta,htheta)

use commondata
use sim_par
use par_size

double precision :: stheta(spx,nz,spy,2), htheta(spx,nz,spy,2)

htheta=dt*stheta

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth_theta(stheta,htheta)

use commondata
use sim_par
use par_size
use sterms

double precision :: stheta(spx,nz,spy,2), htheta(spx,nz,spy,2)

htheta=0.5d0*dt*(3.0d0*stheta-stheta_o)

return
end
