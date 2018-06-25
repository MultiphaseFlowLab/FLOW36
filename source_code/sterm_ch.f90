subroutine sterm_ch(sphi)

use commondata
use par_size
use velocity
use phase_field
use wavenumber
use grid

double precision :: sphi(spx,nz,spy,2),epsnum,mask
double precision :: modnabphi,normflux,funflux(nz)
double precision, allocatable, dimension(:,:,:,:) :: a1,a2,a3
double precision, allocatable, dimension(:,:,:) :: a1f,convf
!USED ONLY WHEN THE CORRECTION IS ENABLED
double precision, allocatable, dimension(:,:,:) :: phinx,phiny,phinz
double precision, allocatable, dimension(:,:,:) :: a4f,a5f,a6f,nabcf
double precision, allocatable, dimension(:,:,:,:) :: nabc,a4,a5,a6


integer :: i,j,k

#define phicorflag phicorcompflag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 0
! STANDARD MODEL

! calculate -(1+s)/Pe*nabla^2 phi
allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

sphi=-(1.0d0+s_coeff)/pe*sphi

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
sphi=sphi-a1

!Calculate nabla^2 \phi^3
a1f=phi**(3.0d0)

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a3,a1,a1f)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 1
! PROFILE-CORRECTED
! calculate -(1+s)/Pe*nabla^2 phi
allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

sphi=-(1.0d0+s_coeff-lamphi)/pe*sphi

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

allocate(phinx(nx,fpz,fpy))
phinx=a1f

! v*(d phi /dy)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

allocate(phiny(nx,fpz,fpy))
phiny=a1f

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

allocate(phinz(nx,fpz,fpy))
phinz=a1f

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
sphi=sphi-a1


!CALCULATE THE LAMBDA TERM (EXPLICIT)
a1f=1.0d0-phi**(2.0d0)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=sqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))


do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a1f(i,k,j)*phinx(i,k,j)
      a5f(i,k,j)=a1f(i,k,j)*phiny(i,k,j)
      a6f(i,k,j)=a1f(i,k,j)*phinz(i,k,j)
    enddo
  enddo
enddo

deallocate(phinx,phiny,phinz)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call phys_to_spectral(a4f,a1,1)
call phys_to_spectral(a5f,a2,1)
call phys_to_spectral(a6f,a3,1)

deallocate(a4f,a5f,a6f)

allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

!X DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a6)

!SUM EVEYTHING TO SPHI
sphi=sphi-(lamphi/(sqrt(2.0d0)*ch*pe))*(a4+a5+a6)

deallocate(a4,a5,a6)

!Calculate nabla^2 \phi^3
a1f=phi**(3.0d0)

call phys_to_spectral(a1f,a1,1)

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a3,a1,a1f)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 2
! Flux-Corrected

allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

allocate(nabc(spx,nz,spy,2))
nabc=ch*ch*sphi
sphi=-(s_coeff-lamphi)/pe*sphi

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinx(nx,fpz,fpy))
phinx=a1f

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)
allocate(phiny(nx,fpz,fpy))
phiny=a1f

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinz(nx,fpz,fpy))
phinz=a1f

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(convf,a1,1)
deallocate(convf)
! sum all the convective terms to sphi
sphi=sphi-a1

! computing 1/Pe nabla^2 phi^3
a1f=phi**(3d0)

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a3)

!COMPUTING THE DIVERGENCE TERM..
nabc=-a1+nabc

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*nabc(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*nabc(i,:,j,1)
  enddo
enddo

do j=1,spy
  do i=1,spx
    a2(i,:,j,1)=-ky(j+cstart(3))*nabc(i,:,j,2)
    a2(i,:,j,2)=ky(j+cstart(3))*nabc(i,:,j,1)
  enddo
enddo

call dz(nabc,a3)

deallocate(nabc)

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a1,a4f,1)
call spectral_to_phys(a2,a5f,1)
call spectral_to_phys(a3,a6f,1)

a1f=1d0-phi**(2d0)

!!NORMALIZATION OF THE PHI GRADIENT
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=sqrt(phinx(i,k,j)**2d0+phiny(i,k,j)**2d0+phinz(i,k,j)**2d0)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      normflux=a4f(i,k,j)*phinx(i,k,j)+a5f(i,k,j)*phiny(i,k,j)+a6f(i,k,j)*phinz(i,k,j)
      a4f(i,k,j)=-lamphi/(sqrt(2d0)*ch)*a1f(i,k,j)*phinx(i,k,j) + normflux*phinx(i,k,j)
      a5f(i,k,j)=-lamphi/(sqrt(2d0)*ch)*a1f(i,k,j)*phiny(i,k,j) + normflux*phiny(i,k,j)
      a6f(i,k,j)=-lamphi/(sqrt(2d0)*ch)*a1f(i,k,j)*phinz(i,k,j) + normflux*phinz(i,k,j)
    enddo
  enddo
enddo

deallocate(a1f,phinx,phiny,phinz)

call phys_to_spectral(a4f,a1,1)
call phys_to_spectral(a5f,a2,1)
call phys_to_spectral(a6f,a3,1)

! COMPUITING THE DIVERGENCE OF THE CORRECTION
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

!X DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a6)

deallocate(a1,a2,a3)

sphi=sphi+(a4+a5+a6)/pe

deallocate(a4,a5,a6)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 3
! PROFILE-CORRECTED DISABLED AT THE WALLS MF!!
! calculate -(1+s)/Pe*nabla^2 phi
allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

allocate(nabc(spx,nz,spy,2))
nabc=sphi
sphi=-(1.0d0+s_coeff)/pe*sphi

!Calculate convective term
allocate(a1f(nx,fpz,fpy))
! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

allocate(phinx(nx,fpz,fpy))
phinx=a1f

! v*(d phi /dy)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

allocate(phiny(nx,fpz,fpy))
phiny=a1f

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

allocate(phinz(nx,fpz,fpy))
phinz=a1f

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
sphi=sphi-a1


!CALCULATE THE LAMBDA TERM (EXPLICIT)
a1f=1.0d0-phi**(2.0d0)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=sqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))


do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a1f(i,k,j)*phinx(i,k,j)
      a5f(i,k,j)=a1f(i,k,j)*phiny(i,k,j)
      a6f(i,k,j)=a1f(i,k,j)*phinz(i,k,j)
    enddo
  enddo
enddo

deallocate(phinx,phiny,phinz)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call phys_to_spectral(a4f,a1,1)
call phys_to_spectral(a5f,a2,1)
call phys_to_spectral(a6f,a3,1)

deallocate(a4f,a5f,a6f)

allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

!X DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a6)

!SUM EVEYTHING TO SPHI
!! DISABLE AT THE WALLS funflux(k)
do k=1,fpz
   funflux(k)=0.5d0*(tanh((abs(z(fstart(2)+k)) - 0.975d0)/0.01d0)+1d0)
!   funflux(k)=0d0
!   if ((fstart(2)+k) .eq. 1)  funflux(k)=1d0
!   if ((fstart(2)+k) .eq. nz) funflux(k)=1d0
!   print*,'RANK',RANK,'funflux',funflux(k)
enddo
!!!
nabc=(lamphi/pe)*nabc-(lamphi/(sqrt(2d0)*ch*pe))*(a4+a5+a6)
deallocate(a4,a5,a6)
allocate(nabcf(nx,fpz,fpy))
call spectral_to_phys(nabc,nabcf,1)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      nabcf(i,k,j)=(1d0-funflux(k))*nabcf(i,k,j)
    enddo
  enddo
enddo
call phys_to_spectral(nabcf,nabc,1)
deallocate(nabcf)

sphi=sphi+nabc
deallocate(nabc)

!Calculate nabla^2 \phi^3
a1f=phi**(3.0d0)

call phys_to_spectral(a1f,a1,1)

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a3,a1,a1f)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 4
! PROFILE CORRECTION kill the gradients
! calculate -(1+s)/Pe*nabla^2 phi
allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

sphi=-(1.0d0+s_coeff-lamphi)/pe*sphi

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

allocate(phinx(nx,fpz,fpy))
phinx=a1f

! v*(d phi /dy)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

allocate(phiny(nx,fpz,fpy))
phiny=a1f

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

allocate(phinz(nx,fpz,fpy))
phinz=a1f

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
sphi=sphi-a1


! calculate nabla^2 phi^3
a1f=phi**(3.0d0)

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))


call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a3)


! calculate profile-corrected flux
a1f=1.0d0-phi**(2.0d0)

epsnum=1.0d0/(50.0d0*Ch)

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=sqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      mask=max(modnabphi-epsnum,0d0)/(modnabphi-epsnum)
      a4f(i,k,j)=-1.0d0/(sqrt(2.0d0)*Ch)*a1f(i,k,j)*phinx(i,k,j)/modnabphi
      a5f(i,k,j)=-1.0d0/(sqrt(2.0d0)*Ch)*a1f(i,k,j)*phiny(i,k,j)/modnabphi
      a6f(i,k,j)=-1.0d0/(sqrt(2.0d0)*Ch)*a1f(i,k,j)*phinz(i,k,j)/modnabphi
      a4f(i,k,j)=a4f(i,k,j)*mask
      a5f(i,k,j)=a5f(i,k,j)*mask
      a6f(i,k,j)=a6f(i,k,j)*mask
    enddo
  enddo
enddo

deallocate(phinx,phiny,phinz,a1f)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call phys_to_spectral(a4f,a1,1)
call phys_to_spectral(a5f,a2,1)
call phys_to_spectral(a6f,a3,1)

deallocate(a4f,a5f,a6f)

allocate(a4(spx,nz,spy,2))

!X DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a2)

!SUM EVEYTHING TO SPHI
sphi=sphi+(lamphi/pe)*(a4+a1+a2)

deallocate(a1,a2,a3,a4)

#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 5
! Flux-Corrected kill the gradients

allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

allocate(nabc(spx,nz,spy,2))
nabc=ch**2.0d0*sphi
sphi=-(1.0d0+s_coeff-lamphi)/pe*sphi

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinx(nx,fpz,fpy))
phinx=a1f

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)
allocate(phiny(nx,fpz,fpy))
phiny=a1f

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinz(nx,fpz,fpy))
phinz=a1f

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(convf,a1,1)
deallocate(convf)
! sum all the convective terms to sphi
sphi=sphi-a1

! computing 1/Pe nabla^2 phi^3
a1f=phi**(3d0)

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a3)

!COMPUTING THE DIVERGENCE TERM..
nabc=-a1+phic+nabc

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*nabc(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*nabc(i,:,j,1)
  enddo
enddo

do j=1,spy
  do i=1,spx
    a2(i,:,j,1)=-ky(j+cstart(3))*nabc(i,:,j,2)
    a2(i,:,j,2)=ky(j+cstart(3))*nabc(i,:,j,1)
  enddo
enddo

call dz(nabc,a3)

deallocate(nabc)

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a1,a4f,1)
call spectral_to_phys(a2,a5f,1)
call spectral_to_phys(a3,a6f,1)

a1f=1d0-phi**(2d0)

!!NORMALIZATION OF THE PHI GRADIENT
epsnum=1.0d0/(50.0d0*Ch)


do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=sqrt(phinx(i,k,j)**2d0+phiny(i,k,j)**2d0+phinz(i,k,j)**2d0)
      mask=max(modnabphi-epsnum,0d0)/(modnabphi-epsnum)
      normflux=(a4f(i,k,j)*phinx(i,k,j)+a5f(i,k,j)*phiny(i,k,j)+a6f(i,k,j)*phinz(i,k,j))/modnabphi
      a4f(i,k,j)=+(-lamphi/(sqrt(2d0)*ch)*a1f(i,k,j)+normflux)*phinx(i,k,j)/modnabphi
      a5f(i,k,j)=+(-lamphi/(sqrt(2d0)*ch)*a1f(i,k,j)+normflux)*phiny(i,k,j)/modnabphi
      a6f(i,k,j)=+(-lamphi/(sqrt(2d0)*ch)*a1f(i,k,j)+normflux)*phinz(i,k,j)/modnabphi
      a4f(i,k,j)=a4f(i,k,j)*mask
      a5f(i,k,j)=a5f(i,k,j)*mask
      a6f(i,k,j)=a6f(i,k,j)*mask
    enddo
  enddo
enddo

deallocate(a1f,phinx,phiny,phinz)

call phys_to_spectral(a4f,a1,1)
call phys_to_spectral(a5f,a2,1)
call phys_to_spectral(a6f,a3,1)

! COMPUITING THE DIVERGENCE OF THE CORRECTION
allocate(a4(spx,nz,spy,2))

!X DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a2)



sphi=sphi+(a4+a1+a2)/pe

deallocate(a1,a2,a3,a4)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

return
end
