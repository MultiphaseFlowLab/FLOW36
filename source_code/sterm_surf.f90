subroutine sterm_surf(spsi)

use commondata
use par_size
use velocity
use phase_field
use surfactant
use wavenumber
use dual_grid

double precision :: spsi(spxpsi,npsiz,spypsi,2)
double precision, dimension(npsix,fpzpsi,fpypsi) :: a1,a2,a3,a4,fspsi
double precision :: b1s(spxpsi,npsiz,spypsi,2)

integer :: i,j,k

call coarse2fine(phic,phic_fg)
call spectral_to_phys_fg(phic_fg,phi_fg,1,1)
call spectral_to_phys_fg(psic_fg,psi_fg,1,1)


do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      a1(i,k,j)=psi_fg(i,k,j)*(1.0d0-psi_fg(i,k,j))*(2.0d0*phi_fg(i,k,j)-2.0d0*phi_fg(i,k,j)**3+phi_fg(i,k,j)/Ex)
      a2(i,k,j)=psi_fg(i,k,j)*(1.0d0-psi_fg(i,k,j))*(2.0d0 - 6.0d0*phi_fg(i,k,j)**2 + 1.0d0/Ex)
      a3(i,k,j)=(1.0d0-2.0d0*psi_fg(i,k,j))        *(2.0d0*phi_fg(i,k,j)-2.0d0*phi_fg(i,k,j)**3+phi_fg(i,k,j)/Ex)
    enddo
  enddo
enddo

a1=a1/pe_psi
a2=a2/pe_psi
a3=a3/pe_psi

! laplacian of phi
call dz_fg(phic_fg,b1s)
call dz_fg(b1s,spsi)

do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      spsi(i,k,j,1)=spsi(i,k,j,1)-k2psi(i+cstartpsi(1),j+cstartpsi(3))*phic_fg(i,k,j,1)
      spsi(i,k,j,2)=spsi(i,k,j,2)-k2psi(i+cstartpsi(1),j+cstartpsi(3))*phic_fg(i,k,j,2)
    enddo
  enddo
enddo

call spectral_to_phys_fg(spsi,fspsi,1,1)

do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      fspsi(i,k,j)=a1(i,k,j)*fspsi(i,k,j)
    enddo
  enddo
enddo


! d phi/dx
b1s=0.0d0
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      b1s(i,k,j,1)=-kxpsi(i+cstartpsi(1))*phic_fg(i,k,j,2)
      b1s(i,k,j,2)=kxpsi(i+cstartpsi(1))*phic_fg(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys_fg(b1s,a1,1,1)   ! d phi/dx

b1s=0.0d0
a4=0.0d0
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      b1s(i,k,j,1)=-kxpsi(i+cstartpsi(1))*psic_fg(i,k,j,2)
      b1s(i,k,j,2)=kxpsi(i+cstartpsi(1))*psic_fg(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys_fg(b1s,a4,1,1)   ! d psi/dx

do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      fspsi(i,k,j)=fspsi(i,k,j)+a2(i,k,j)*a1(i,k,j)**2 + a3(i,k,j)*a1(i,k,j)*a4(i,k,j) - u_fg(i,k,j)*a4(i,k,j)
    enddo
  enddo
enddo


! d phi/dy
b1s=0.0d0
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      b1s(i,k,j,1)=-kypsi(j+cstartpsi(3))*phic_fg(i,k,j,2)
      b1s(i,k,j,2)=kypsi(j+cstartpsi(3))*phic_fg(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys_fg(b1s,a1,1,1)   ! d phi/dy

b1s=0.0d0
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      b1s(i,k,j,1)=-kypsi(j+cstartpsi(3))*psic_fg(i,k,j,2)
      b1s(i,k,j,2)=kypsi(j+cstartpsi(3))*psic_fg(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys_fg(b1s,a4,1,1)   ! d psi/dy

do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      fspsi(i,k,j)=fspsi(i,k,j)+a2(i,k,j)*a1(i,k,j)**2 + a3(i,k,j)*a1(i,k,j)*a4(i,k,j) - v_fg(i,k,j)*a4(i,k,j)
    enddo
  enddo
enddo


! d phi/dz
call dz_fg(phic_fg,b1s)

call spectral_to_phys_fg(b1s,a1,1,1)   ! d phi/dz

call dz_fg(psic_fg,b1s)

call spectral_to_phys_fg(b1s,a4,1,1)   ! d psi/dz

do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      fspsi(i,k,j)=fspsi(i,k,j)+a2(i,k,j)*a1(i,k,j)**2 + a3(i,k,j)*a1(i,k,j)*a4(i,k,j) - w_fg(i,k,j)*a4(i,k,j)
    enddo
  enddo
enddo
! transform non-linear term to spectral space
call phys_to_spectral_fg(fspsi,spsi,1,1)

return
end
