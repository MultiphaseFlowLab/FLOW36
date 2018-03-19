subroutine chop_modes(varc)

use commondata
use par_size

double precision :: varc(spx,nz,spy,2)

integer :: limx,llimy,ulimy

! chop high Fourier x modes
limx=floor(2.0/3.0*real(nx/2+1))+1
if(cstart(1)+spx.ge.limx)then
!  if(cstart(1).ge.limx)then
!    varc(1:spx,1:nz,1:spy,1)=0.0d0
!    varc(1:spx,1:nz,1:spy,2)=0.0d0
!  else
!    varc(limx-cstart(1):spx,1:nz,1:spy,1)=0.0d0
!    varc(limx-cstart(1):spx,1:nz,1:spy,2)=0.0d0
!  endif
  varc(max(limx-cstart(1),1):spx,1:nz,1:spy,1)=0.0d0
  varc(max(limx-cstart(1),1):spx,1:nz,1:spy,2)=0.0d0
endif


! chop high Fourier y modes
llimy=floor(2.0/3.0*real(ny/2+1))
ulimy=ny-floor(2.0/3.0*real(ny/2))
if((cstart(3)+spy.ge.llimy).and.(cstart(3).le.ulimy))then
  varc(1:spx,1:nz,max(llimy-cstart(3),1):min(ulimy-cstart(3),spy),1)=0.0d0
  varc(1:spx,1:nz,max(llimy-cstart(3),1):min(ulimy-cstart(3),spy),2)=0.0d0
endif


! chop high Chebyshev modes
varc(1:spx,floor(2.0*real(nz)/3.0)+1:nz,1:spy,1)=0.0d0
varc(1:spx,floor(2.0*real(nz)/3.0)+1:nz,1:spy,2)=0.0d0

return
end
