subroutine chop_modes(varc)

use commondata
use par_size

double precision :: varc(spx,nz,spy,2)

integer :: limx,llimy,ulimy

! chop high Fourier x modes
!$acc kernels
limx=floor(2.0/3.0*real(nx/2+1))+1
if(cstart(1)+spx.ge.limx)then
  varc(max(limx-cstart(1),1):spx,1:nz,1:spy,1)=0.0d0
  varc(max(limx-cstart(1),1):spx,1:nz,1:spy,2)=0.0d0
endif
!$acc end kernels

! chop high Fourier y modes
!$acc kernels
llimy=floor(2.0/3.0*real(ny/2+1))+1
ulimy=ny-floor(2.0/3.0*real(ny/2))
if((cstart(3)+spy.ge.llimy).and.(cstart(3).le.ulimy))then
  varc(1:spx,1:nz,max(llimy-cstart(3),1):min(ulimy-cstart(3),spy),1)=0.0d0
  varc(1:spx,1:nz,max(llimy-cstart(3),1):min(ulimy-cstart(3),spy),2)=0.0d0
endif
!$acc end kernels

! chop high Chebyshev modes
!$acc kernels
varc(1:spx,floor(2.0*real(nz)/3.0)+1:nz,1:spy,1)=0.0d0
varc(1:spx,floor(2.0*real(nz)/3.0)+1:nz,1:spy,2)=0.0d0
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chop_modes_fg(varc)

use commondata
use par_size
use dual_grid

double precision :: varc(spxpsi,npsiz,spypsi,2)

integer :: limx,llimy,ulimy

! chop high Fourier x modes
!$acc kernels
limx=floor(2.0/3.0*real(npsix/2+1))+1
if(cstartpsi(1)+spxpsi.ge.limx)then
  varc(max(limx-cstartpsi(1),1):spxpsi,1:npsiz,1:spypsi,1)=0.0d0
  varc(max(limx-cstartpsi(1),1):spxpsi,1:npsiz,1:spypsi,2)=0.0d0
endif
!$acc end kernels


! chop high Fourier y modes
!$acc kernels
llimy=floor(2.0/3.0*real(npsiy/2+1))+1
ulimy=npsiy-floor(2.0/3.0*real(npsiy/2))
if((cstartpsi(3)+spypsi.ge.llimy).and.(cstartpsi(3).le.ulimy))then
  varc(1:spxpsi,1:npsiz,max(llimy-cstartpsi(3),1):min(ulimy-cstartpsi(3),spypsi),1)=0.0d0
  varc(1:spxpsi,1:npsiz,max(llimy-cstartpsi(3),1):min(ulimy-cstartpsi(3),spypsi),2)=0.0d0
endif
!$acc end kernels


! chop high Chebyshev modes
!$acc kernels
varc(1:spxpsi,floor(2.0*real(npsiz)/3.0)+1:npsiz,1:spypsi,1)=0.0d0
varc(1:spxpsi,floor(2.0*real(npsiz)/3.0)+1:npsiz,1:spypsi,2)=0.0d0
!$acc end kernels

return
end
