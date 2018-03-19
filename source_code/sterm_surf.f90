subroutine sterm_surf(spsi)

use commondata
use par_size
use velocity
use phase_field
use surfactant
use wavenumber

double precision :: spsi(spx,nz,spy,2)
double precision :: a1(nx,fpz,fpy),a2(nx,fpz,fpy),a3(nx,fpz,fpy),a4(nx,fpz,fpy),fspsi(nx,fpz,fpy)
double precision :: b1s(spx,nz,spy,2)

integer :: i,j,k

call spectral_to_phys(phic,phi,1)
call spectral_to_phys(psic,psi,1)


do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a1(i,k,j)=psi(i,k,j)*(1.0d0-psi(i,k,j))*(2.0d0*phi(i,k,j)-2.0d0*phi(i,k,j)**3+phi(i,k,j)/Ex)
      a2(i,k,j)=psi(i,k,j)*(1.0d0-psi(i,k,j))*(2.0d0 - 6.0d0*phi(i,k,j)**2 + 1.0d0/Ex)
      a3(i,k,j)=(1.0d0-2.0d0*psi(i,k,j))     *(2.0d0*phi(i,k,j)-2.0d0*phi(i,k,j)**3+phi(i,k,j)/Ex)
    enddo
  enddo
enddo

a1=a1/pe_psi
a2=a2/pe_psi
a3=a3/pe_psi

! laplacian of phi
call dz(phic,b1s)
call dz(b1s,spsi)

do j=1,spy
  do k=1,nz
    do i=1,spx
      spsi(i,k,j,1)=spsi(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,k,j,1)
      spsi(i,k,j,2)=spsi(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,k,j,2)
    enddo
  enddo
enddo

call spectral_to_phys(spsi,fspsi,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      fspsi(i,k,j)=a1(i,k,j)*fspsi(i,k,j)
    enddo
  enddo
enddo


! d phi/dx
b1s=0.0d0
do j=1,spy
  do k=1,nz
    do i=1,spx
      b1s(i,k,j,1)=-kx(i+cstart(1))*phic(i,k,j,2)
      b1s(i,k,j,2)=kx(i+cstart(1))*phic(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(b1s,a1,1)   ! d phi/dx

call spectral_to_phys(uc,u,1)

b1s=0.0d0
a4=0.0d0
do j=1,spy
  do k=1,nz
    do i=1,spx
      b1s(i,k,j,1)=-kx(i+cstart(1))*psic(i,k,j,2)
      b1s(i,k,j,2)=kx(i+cstart(1))*psic(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(b1s,a4,1)   ! d psi/dx

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      fspsi(i,k,j)=fspsi(i,k,j)+a2(i,k,j)*a1(i,k,j)**2 + a3(i,k,j)*a1(i,k,j)*a4(i,k,j) - u(i,k,j)*a4(i,k,j)
    enddo
  enddo
enddo


! d phi/dy
b1s=0.0d0
do j=1,spy
  do k=1,nz
    do i=1,spx
      b1s(i,k,j,1)=-ky(j+cstart(3))*phic(i,k,j,2)
      b1s(i,k,j,2)=ky(j+cstart(3))*phic(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(b1s,a1,1)
   ! d phi/dy
call spectral_to_phys(vc,v,1)

b1s=0.0d0
do j=1,spy
  do k=1,nz
    do i=1,spx
      b1s(i,k,j,1)=-ky(j+cstart(3))*psic(i,k,j,2)
      b1s(i,k,j,2)=ky(j+cstart(3))*psic(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys(b1s,a4,1)   ! d psi/dy

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      fspsi(i,k,j)=fspsi(i,k,j)+a2(i,k,j)*a1(i,k,j)**2 + a3(i,k,j)*a1(i,k,j)*a4(i,k,j) - v(i,k,j)*a4(i,k,j)
    enddo
  enddo
enddo


! d phi/dz
call dz(phic,b1s)

call spectral_to_phys(b1s,a1,1)   ! d phi/dz

call spectral_to_phys(wc,w,1)

call dz(psic,b1s)

call spectral_to_phys(b1s,a4,1)   ! d psi/dz

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      fspsi(i,k,j)=fspsi(i,k,j)+a2(i,k,j)*a1(i,k,j)**2 + a3(i,k,j)*a1(i,k,j)*a4(i,k,j) - w(i,k,j)*a4(i,k,j)
    enddo
  enddo
enddo
! transform non-linear term to spectral space
call phys_to_spectral(fspsi,spsi,1)

return
end
