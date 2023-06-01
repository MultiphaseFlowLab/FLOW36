subroutine sterm_ch(sphi)

use commondata
use par_size
use velocity
use phase_field
use wavenumber
use grid

double precision :: sphi(spx,nz,spy,2),epsnum,mask,epsnum6
double precision :: modnabphi,normflux,funflux(nz)
double precision :: mink,maxk
double precision, allocatable, dimension(:,:,:,:) :: a1,a2,a3
double precision, allocatable, dimension(:,:,:) :: a1f,a2f,a3f,convf
!use only when phicor_flag !=0 
double precision, allocatable, dimension(:,:,:) :: phinx,phiny,phinz,phinmod
double precision, allocatable, dimension(:,:,:) :: a4f,a5f,a6f,nabcf
double precision, allocatable, dimension(:,:,:,:) :: nabc,a4,a5,a6



integer :: i,j,k

#define phicorflag phicorcompflag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sphi=0.0d0

#if phicorflag == 0
! Standard MODEL

! calculate -(1+s)/Pe*nabla^2 phi
allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

!$acc kernels
sphi=-(1.0d0+s_coeff)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(2)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
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

!$acc parallel loop collapse(3)
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
!$acc kernels
sphi=sphi-a1

!Calculate nabla^2 \phi^3
a1f=phi**(3.0d0)
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

!$acc parallel loop collapse(3)
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

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

!$acc kernels
sphi=-(1.0d0+s_coeff-lamphi)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

allocate(phiny(nx,fpz,fpy))

!$acc kernels
phiny=a1f
!$acc end kernels

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

allocate(phinz(nx,fpz,fpy))
!$acc kernels
phinz=a1f
!$acc end kernels

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
!$acc kernels
sphi=sphi-a1
!$acc end kernels

!CALCULATE THE LAMBDA TERM (EXPLICIT)
!$acc kernels
a1f=1.0d0-phi**(2.0d0)
!$acc end kernels

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=dsqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

!$acc parallel loop collapse(3)
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
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a6)

!SUM EVEYTHING TO SPHI
!$acc kernels
sphi=sphi-(lamphi/(dsqrt(2.0d0)*ch*pe))*(a4+a5+a6)
!$acc end kernels

deallocate(a4,a5,a6)

!Calculate nabla^2 \phi^3
!$acc kernels
a1f=phi**(3.0d0)
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

!$acc parallel loop collapse(3)
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

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

allocate(nabc(spx,nz,spy,2))
!$acc kernels
nabc=ch*ch*sphi
sphi=-(s_coeff-lamphi)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)
allocate(phiny(nx,fpz,fpy))
!$acc kernels
phiny=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
!$acc kernels
phinz=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
!$acc kernels
sphi=sphi-a1
!$acc end kernels

! computing 1/Pe nabla^2 phi^3
!$acc kernels
a1f=phi**(3d0)
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

!$acc parallel loop collapse(3)
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
!$acc kernels
nabc=-a1+nabc
!$acc end kernels

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*nabc(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*nabc(i,:,j,1)
  enddo
enddo

!$acc parallel loop collapse(2)
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

!$acc kernels
a1f=1d0-phi**(2d0)
!$acc end kernels

!!NORMALIZATION OF THE PHI GRADIENT
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=dsqrt(phinx(i,k,j)**2d0+phiny(i,k,j)**2d0+phinz(i,k,j)**2d0)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      normflux=a4f(i,k,j)*phinx(i,k,j)+a5f(i,k,j)*phiny(i,k,j)+a6f(i,k,j)*phinz(i,k,j)
      a4f(i,k,j)=-lamphi/(dsqrt(2d0)*ch)*a1f(i,k,j)*phinx(i,k,j) + normflux*phinx(i,k,j)
      a5f(i,k,j)=-lamphi/(dsqrt(2d0)*ch)*a1f(i,k,j)*phiny(i,k,j) + normflux*phiny(i,k,j)
      a6f(i,k,j)=-lamphi/(dsqrt(2d0)*ch)*a1f(i,k,j)*phinz(i,k,j) + normflux*phinz(i,k,j)
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
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a6)

deallocate(a1,a2,a3)

!$acc kernels
sphi=sphi+(a4+a5+a6)/pe
!$acc end kernels

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

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

allocate(nabc(spx,nz,spy,2))
!$acc kernels
nabc=sphi
sphi=-(1.0d0+s_coeff)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))
! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

allocate(phiny(nx,fpz,fpy))
!$acc kernels
phiny=a1f
!$acc end kernels

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

allocate(phinz(nx,fpz,fpy))
!$acc kernels
phinz=a1f
!$acc end kernels

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
!$acc kernels
sphi=sphi-a1
!$acc end kernels


!CALCULATE THE LAMBDA TERM (EXPLICIT)
!$acc kernels
a1f=1.0d0-phi**(2.0d0)
!$acc end kernels
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=dsqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

!$acc parallel loop collapse(3)
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
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
!$acc parallel loop collapse(2)
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
!$acc parallel loop
do k=1,fpz
   funflux(k)=0.5d0*(tanh((abs(z(fstart(2)+k)) - 0.975d0)/0.01d0)+1d0)
!   funflux(k)=0d0
!   if ((fstart(2)+k) .eq. 1)  funflux(k)=1d0
!   if ((fstart(2)+k) .eq. nz) funflux(k)=1d0
!   print*,'RANK',RANK,'funflux',funflux(k)
enddo
!!!
!$acc kernels
nabc=(lamphi/pe)*nabc-(lamphi/(dsqrt(2d0)*ch*pe))*(a4+a5+a6)
!$acc end kernels
deallocate(a4,a5,a6)
allocate(nabcf(nx,fpz,fpy))
call spectral_to_phys(nabc,nabcf,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      nabcf(i,k,j)=(1d0-funflux(k))*nabcf(i,k,j)
    enddo
  enddo
enddo
call phys_to_spectral(nabcf,nabc,1)
deallocate(nabcf)

!$acc kernels
sphi=sphi+nabc
!$acc end kernels

deallocate(nabc)

!Calculate nabla^2 \phi^3
!$acc kernels
a1f=phi**(3.0d0)
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

!$acc parallel loop collapse(2)
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

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

!$acc kernels
sphi=-(1.0d0+s_coeff-lamphi)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*v(i,k,j)
    enddo
  enddo
enddo

allocate(phiny(nx,fpz,fpy))
!$acc kernels
phiny=a1f
!$acc end kernels

! w*(d phi /dz)
call dz(phic,a1)

call spectral_to_phys(wc,w,1)
call spectral_to_phys(a1,a1f,1)

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=convf(i,k,j)+a1f(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo

allocate(phinz(nx,fpz,fpy))
!$acc kernels
phinz=a1f
!$acc end kernels

call phys_to_spectral(convf,a1,1)

deallocate(convf)
! sum all the convective terms to sphi
!$acc kernels
sphi=sphi-a1
!$acc end kernels

! calculate nabla^2 phi^3
!$acc kernels
a1f=phi**(3.0d0)
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))


call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

!$acc parallel loop collapse(3)
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
!$acc kernels
a1f=1.0d0-phi**(2.0d0)
!$acc end kernels

!$acc kernels
epsnum=1.0d0/(50.0d0*Ch)
!$acc end kernels

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=dsqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      mask=max(modnabphi-epsnum,0d0)/(modnabphi-epsnum)
      a4f(i,k,j)=-1.0d0/(dsqrt(2.0d0)*Ch)*a1f(i,k,j)*phinx(i,k,j)/modnabphi
      a5f(i,k,j)=-1.0d0/(dsqrt(2.0d0)*Ch)*a1f(i,k,j)*phiny(i,k,j)/modnabphi
      a6f(i,k,j)=-1.0d0/(dsqrt(2.0d0)*Ch)*a1f(i,k,j)*phinz(i,k,j)/modnabphi
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
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a2)

!SUM EVEYTHING TO SPHI
!$acc kernels
sphi=sphi+(lamphi/pe)*(a4+a1+a2)
!$acc end kernels

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

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

allocate(nabc(spx,nz,spy,2))
!$acc kernels
nabc=ch**2.0d0*sphi
sphi=-(1.0d0+s_coeff-lamphi)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

allocate(convf(nx,fpz,fpy))

do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)
allocate(phiny(nx,fpz,fpy))
!$acc kernels
phiny=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
!$acc kernels
phinz=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
!$acc kernels
sphi=sphi-a1
!$acc end kernels

! computing 1/Pe nabla^2 phi^3
!$acc kernels
a1f=phi**(3d0)
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

call dz(a1,a2)
call dz(a2,a3)

deallocate(a2)

!$acc parallel loop collapse(3)
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
!$acc kernels
nabc=-a1+phic+nabc
!$acc end kernels

allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*nabc(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*nabc(i,:,j,1)
  enddo
enddo

!$acc parallel loop collapse(2)
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

!$acc kernels
a1f=1d0-phi**(2d0)
!$acc end kernels

!!NORMALIZATION OF THE PHI GRADIENT
epsnum=1.0d0/(50.0d0*Ch)


!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=dsqrt(phinx(i,k,j)**2d0+phiny(i,k,j)**2d0+phinz(i,k,j)**2d0)
      mask=max(modnabphi-epsnum,0d0)/(modnabphi-epsnum)
      normflux=(a4f(i,k,j)*phinx(i,k,j)+a5f(i,k,j)*phiny(i,k,j)+a6f(i,k,j)*phinz(i,k,j))/modnabphi
      a4f(i,k,j)=+(-lamphi/(dsqrt(2d0)*ch)*a1f(i,k,j)+normflux)*phinx(i,k,j)/modnabphi
      a5f(i,k,j)=+(-lamphi/(dsqrt(2d0)*ch)*a1f(i,k,j)+normflux)*phiny(i,k,j)/modnabphi
      a6f(i,k,j)=+(-lamphi/(dsqrt(2d0)*ch)*a1f(i,k,j)+normflux)*phinz(i,k,j)/modnabphi
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
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!Y DERIVATIVE OF CORRECTION
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!Z DERIVATIVE OF CORRECTION
call dz(a3,a2)


!$acc end kernels
sphi=sphi+(a4+a1+a2)/pe
!$acc end kernels

deallocate(a1,a2,a3,a4)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if phicorflag == 6
! Kwakkel model (A redefined energy functional to prevent mass loss in phase-field methods)

! calculate -(1+s)/Pe*nabla^2 phi
allocate(a1(spx,nz,spy,2))

call dz(phic,a1)
call dz(a1,sphi)

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sphi(i,:,j,1)=sphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    sphi(i,:,j,2)=sphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo

!$acc kernels
sphi=-(1.0d0+s_coeff)/pe*sphi
!$acc end kernels

!Calculate convective term
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)
allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(2)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)
allocate(phiny(nx,fpz,fpy))
!$acc kernels
phiny=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
!$acc kernels
phinz=a1f
!$acc end kernels


!$acc parallel loop collapse(3)
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
!$acc kernels
sphi=sphi-a1
!$acc end kernels


!Calculate nabla^2 ( \phi^3 - Ch^2 div*(nablaphi/modnabphi * modnabphi)

allocate(phinmod(nx,fpz,fpy))
allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))
allocate(a4(spx,nz,spy,2))
allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

! avoid NaN on 1/phinmod
epsnum=0.1d0 

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      if (phinx(i,k,j) .le. epsnum) phinx(i,k,j)=0.0d0
      if (phiny(i,k,j) .le. epsnum) phiny(i,k,j)=0.0d0
      if (phinz(i,k,j) .le. epsnum) phinz(i,k,j)=0.0d0
      phinmod(i,k,j)=dsqrt(phinx(i,k,j)**2d0+phiny(i,k,j)**2d0+phinz(i,k,j)**2d0)
      if (phinmod(i,k,j) .ge. epsnum) then
      a4f(i,k,j)=phinx(i,k,j)/(phinmod(i,k,j))
      a5f(i,k,j)=phiny(i,k,j)/(phinmod(i,k,j))
      a6f(i,k,j)=phinz(i,k,j)/(phinmod(i,k,j))
      else
        a4f(i,k,j)=0.0d0
        a5f(i,k,j)=0.0d0
        a6f(i,k,j)=0.0d0
      endif
      write(*,*) "phinx", a4f(i,k,j)
    enddo
  enddo
enddo

deallocate(phinx,phiny,phinz)

call phys_to_spectral(a4f,a1,1)
call phys_to_spectral(a5f,a2,1)
call phys_to_spectral(a6f,a3,1)

! computing the divergence of the normal 
deallocate(a4f,a5f,a6f)

!x derivative 
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!y derivative 
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!z derivative 
call dz(a3,a2)

! curvature vector
!$acc kernels
a1=a1+a2+a4
!$acc end kernels

! curvature in physical
call spectral_to_phys(a1,a1f,1)

!! Check min and max of curvature
maxk=0.0d0
mink=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      if (abs(a1f(i,k,j)) .ge. 400.d0) a1f(i,k,j)=0.0d0
    enddo
  enddo
enddo

!! Check min and max of curvature
maxk=0.0d0
mink=0.0d0
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      maxk=max(a1f(i,k,j),maxk)
      mink=min(a1f(i,k,j),mink)
    enddo
  enddo
enddo

write(*,*) "max k + rank", maxk, rank
write(*,*) "min k + rank", mink, rank

!do k=1,fpz
!  funflux(k)=0.5d0*(tanh((abs(z(fstart(2)+k)) - 0.975d0)/0.01d0)+1d0)
!   funflux(k)=0d0
!   if ((fstart(2)+k) .eq. 1)  funflux(k)=1d0
!   if ((fstart(2)+k) .eq. nz) funflux(k)=1d0
!   print*,'RANK',RANK,'funflux',funflux(k)
!enddo

! assembly Ch^2*curvature*phinmod
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
!    if (((fstart(2)+k) .le. 10) .or. (fstart(2)+k .ge. nz-9)) then 
     !mask=1.0d0
     !if (abs(phi(i,k,j)) .ge. 0.97d0) mask=0.0d0
     a1f(i,k,j)=phi(i,k,j)**(3.0d0) - ch*ch*a1f(i,k,j)*phinmod(i,k,j) 
    enddo
  enddo
enddo

! assembly in spectral
call phys_to_spectral(a1f,a1,1)

call dz(a1,a2)
call dz(a2,a3)

deallocate(a1f,a2,a4,phinmod)

!$acc parallel loop collapse(3)
do j=1,spy
  do k=1,nz
    do i=1,spx
      sphi(i,k,j,1)=sphi(i,k,j,1)+(a3(i,k,j,1)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,1))/pe
      sphi(i,k,j,2)=sphi(i,k,j,2)+(a3(i,k,j,2)-k2(i+cstart(1),j+cstart(3))*a1(i,k,j,2))/pe
    enddo
  enddo
enddo

deallocate(a1,a3)

#endif


#if phicorflag == 7
! Allen-Cahn, 2nd order phase-field method, conservatvie version of Mirjalili
! Used in combination with calculate_phi_ac(hphi), see calculate_var.f90

allocate(a1(spx,nz,spy,2))
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)

allocate(phinx(nx,fpz,fpy))
!$acc kernels
phinx=a1f
!$acc end kernels

allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(2)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)

allocate(phiny(nx,fpz,fpy))
!$acc kernels
phiny=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
!$acc kernels
phinz=a1f
!$acc end kernels

!$acc parallel loop collapse(3)
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
sphi= -a1

! compute the corrective term
!$acc kernels
a1f=1.0d0-phi**(2.0d0)
!$acc end kernels

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=dsqrt(phinx(i,k,j)**2+phiny(i,k,j)**2+phinz(i,k,j)**2)
      phinx(i,k,j)=phinx(i,k,j)/modnabphi
      phiny(i,k,j)=phiny(i,k,j)/modnabphi
      phinz(i,k,j)=phinz(i,k,j)/modnabphi
    enddo
  enddo
enddo

allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

!$acc parallel loop collapse(3)
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

!x derivative of correction
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a4(i,:,j,2)=kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo
!y derivative of correction
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a2(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a2(i,:,j,1)
  enddo
enddo
!z derivative of correction
call dz(a3,a6)

deallocate(a1,a1f,a3)

!sum everything to sphi 
!$acc kernels
sphi=sphi - 1.0d0/pe*(a4+a5+a6)
!$acc end kernels

deallocate(a4,a5,a6)

#endif




#if phicorflag == 8
! Allen-Cahn, 2nd order phase-field method, conservatvie version of Jain, improved versio of MIrjalili.
! Used in combination with calculate_phi_ac(hphi), see calculate_var.f90
! Article: Accurate conservative phase-field method for simulation of two-phase flows

allocate(a1(spx,nz,spy,2))
allocate(a1f(nx,fpz,fpy))

! u*(d phi /dx)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(uc,u,1)
call spectral_to_phys(a1,a1f,1)


allocate(convf(nx,fpz,fpy))

!$acc parallel loop collapse(2)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      convf(i,k,j)=a1f(i,k,j)*u(i,k,j)
    enddo
  enddo
enddo

! v*(d phi /dy)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    a1(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

call spectral_to_phys(vc,v,1)
call spectral_to_phys(a1,a1f,1)


!$acc parallel loop collapse(3)
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


!$acc parallel loop collapse(3)
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
sphi= -a1

!Compute the sharpening term
epsnum=1.e-4

! a1/a1f is the auxiliary function (psi in Jain paper)
!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a1f(i,k,j)=ch*log((phi(i,k,j)+epsnum)/(1.0d0-phi(i,k,j)+epsnum))
    enddo
  enddo
enddo
!$acc end kernels

call phys_to_spectral(a1f,a1,1)

! to save the gradient of a1f (psi in the paper of Jain)
allocate(a2(spx,nz,spy,2))
allocate(a3(spx,nz,spy,2))
allocate(a4(spx,nz,spy,2))

! X-derivative in spectral of the auxiliary function
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a2(i,:,j,1)=-kx(i+cstart(1))*a1(i,:,j,2)
    a2(i,:,j,2)= kx(i+cstart(1))*a1(i,:,j,1)
  enddo
enddo

! Y-derivative in spectral of the auxiliary function
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a3(i,:,j,1)=-ky(j+cstart(3))*a1(i,:,j,2)
    a3(i,:,j,2)= ky(j+cstart(3))*a1(i,:,j,1)
  enddo
enddo

!Z-derivative in spectral of the auxiliary function
call dz(a1,a4)

!X,Y,Z-derivative in physical of the auxiliar function
allocate(a2f(nx,fpz,fpy))
allocate(a3f(nx,fpz,fpy))
allocate(a4f(nx,fpz,fpy))

!X,Y,Z-derivative in spectral of the auxiliar function
call spectral_to_phys(a2,a2f,1)
call spectral_to_phys(a3,a3f,1)
call spectral_to_phys(a4,a4f,1)


! Normal in physical space
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      modnabphi=max(epsnum,dsqrt(a2f(i,k,j)**2+a3f(i,k,j)**2+a4f(i,k,j)**2))
      a2f(i,k,j)=a2f(i,k,j)/modnabphi
      a3f(i,k,j)=a3f(i,k,j)/modnabphi
      a4f(i,k,j)=a4f(i,k,j)/modnabphi
    enddo
  enddo
enddo



!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a2f(i,k,j)= -0.25d0*(1.0d0 - dtanh(a1f(i,k,j)/2/ch)**2)*a2f(i,k,j)
      a3f(i,k,j)= -0.25d0*(1.0d0 - dtanh(a1f(i,k,j)/2/ch)**2)*a3f(i,k,j)
      a4f(i,k,j)= -0.25d0*(1.0d0 - dtanh(a1f(i,k,j)/2/ch)**2)*a4f(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(a2f,a2,1)
call phys_to_spectral(a3f,a3,1)
call phys_to_spectral(a4f,a4,1)

deallocate(a2f,a3f,a4f)
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

!x derivative of correction
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a1(i,:,j,1)=-kx(i+cstart(1))*a2(i,:,j,2)
    a1(i,:,j,2)=kx(i+cstart(1))*a2(i,:,j,1)
  enddo
enddo
!y derivative of correction
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    a5(i,:,j,1)=-ky(j+cstart(3))*a3(i,:,j,2)
    a5(i,:,j,2)=ky(j+cstart(3))*a3(i,:,j,1)
  enddo
enddo
!z derivative of correction
call dz(a4,a6)

deallocate(a2,a4,a1f,a3)

!sum everything to sphi 
!$acc kernels
sphi=sphi + 1.0d0/pe*(a1+a5+a6)
!$acc end kernels

deallocate(a1,a5,a6)



#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

return
end
