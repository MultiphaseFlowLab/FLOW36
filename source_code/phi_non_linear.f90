subroutine phi_non_linear(s1,s2,s3)

use commondata
use par_size
use phase_field
use wavenumber
use velocity
use sim_par
use velocity_old
use surfactant
use dual_grid


double precision, dimension(spx,nz,spy,2) :: s1,s2,s3
double precision, allocatable, dimension(:,:,:,:) :: gradphix,gradphiy,gradphiz
double precision, allocatable, dimension(:,:,:,:) :: a4,a5,a6,a7,a8,a9
double precision, allocatable, dimension(:,:,:) :: fgradphix,fgradphiy,fgradphiz
double precision, allocatable, dimension(:,:,:) :: a4f,a5f,a6f,a7f,a8f,a9f,a10f,a11f,a12f,viscnon
double precision, allocatable, dimension(:,:,:) :: sigma
double precision :: phif,gammadot

double precision :: modg,threshold,gthreshold

integer :: i,j,k

#define psiflag psicompflag
#define match_visc matched_viscosity
#define match_dens matched_density
#define b_type buoyancytype
#define bodyflag bodycompflag
#define sgradpflag sgradpcompflag
#define eleflag elecompflag
#define non_new_flag non_newtonian


! phi on coarse grid, physical space needed after surface force calculation
call spectral_to_phys(phic,phi,1)

! calculate surface force
call coarse2fine(phic,phic_fg)
call spectral_to_phys_fg(phic_fg,phi_fg,1)

allocate(gradphix(spxpsi,npsiz,spypsi,2))
allocate(gradphiy(spxpsi,npsiz,spypsi,2))
allocate(gradphiz(spxpsi,npsiz,spypsi,2))
allocate(fgradphix(npsix,fpzpsi,fpypsi))
allocate(fgradphiy(npsix,fpzpsi,fpypsi))
allocate(fgradphiz(npsix,fpzpsi,fpypsi))


! calculate gradient of phi
! x and y component
!$acc parallel loop collapse(3)
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      gradphix(i,k,j,1)=-kxpsi(i+cstartpsi(1))*phic_fg(i,k,j,2)
      gradphix(i,k,j,2)=+kxpsi(i+cstartpsi(1))*phic_fg(i,k,j,1)
      gradphiy(i,k,j,1)=-kypsi(j+cstartpsi(3))*phic_fg(i,k,j,2)
      gradphiy(i,k,j,2)=+kypsi(j+cstartpsi(3))*phic_fg(i,k,j,1)
    enddo
  enddo
enddo

call spectral_to_phys_fg(gradphix,fgradphix,1)
call spectral_to_phys_fg(gradphiy,fgradphiy,1)

! z component
call dz_fg(phic_fg,gradphiz)

call spectral_to_phys_fg(gradphiz,fgradphiz,1)

! phi gradient in physical space in fgradphix,fgradphiy,fgradphiz

allocate(sigma(npsix,fpzpsi,fpypsi))

#if psiflag == 0
! constant fixed surface tension, no surfactant
!$acc kernels
sigma=1.0d0
!$acc end kernels
#elif psiflag == 1
! Langmuir EOS for surface tension
call spectral_to_phys_fg(psic_fg,psi_fg,1)
!$acc kernels
sigma=1.0d0+el*log(1.0d0-psi_fg)
! check if sigma<0.5: avoid unphysical value of surface tension (easier way)
sigma=max(sigma,0.5d0)
!$acc end kernels
#endif


allocate(a4f(npsix,fpzpsi,fpypsi))
allocate(a5f(npsix,fpzpsi,fpypsi))
allocate(a6f(npsix,fpzpsi,fpypsi))
allocate(a4(spxpsi,npsiz,spypsi,2))
allocate(a5(spxpsi,npsiz,spypsi,2))
allocate(a6(spxpsi,npsiz,spypsi,2))


! assemble x component of surface tension force
!$acc parallel loop collapse(3)
do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      a4f(i,k,j)=sigma(i,k,j)*((fgradphiy(i,k,j))**2+(fgradphiz(i,k,j))**2)
      a5f(i,k,j)=-sigma(i,k,j)*fgradphix(i,k,j)*fgradphiy(i,k,j)
      a6f(i,k,j)=-sigma(i,k,j)*fgradphix(i,k,j)*fgradphiz(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral_fg(a4f,gradphix,1)
call phys_to_spectral_fg(a5f,gradphiy,1)
call phys_to_spectral_fg(a6f,gradphiz,1)

call dz_fg(gradphiz,a4)

!$acc parallel loop collapse(3)
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      a4(i,k,j,1)=a4(i,k,j,1)-kxpsi(i+cstartpsi(1))*gradphix(i,k,j,2) &
                           & -kypsi(j+cstartpsi(3))*gradphiy(i,k,j,2)
      a4(i,k,j,2)=a4(i,k,j,2)+kxpsi(i+cstartpsi(1))*gradphix(i,k,j,1) &
                           & +kypsi(j+cstartpsi(3))*gradphiy(i,k,j,1)
    enddo
  enddo
enddo

! assemble y component of surface tension force
!$acc parallel loop collapse(3)
do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      a4f(i,k,j)=-sigma(i,k,j)*fgradphix(i,k,j)*fgradphiy(i,k,j)
      a5f(i,k,j)=sigma(i,k,j)*((fgradphix(i,k,j))**2+(fgradphiz(i,k,j))**2)
      a6f(i,k,j)=-sigma(i,k,j)*fgradphiy(i,k,j)*fgradphiz(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral_fg(a4f,gradphix,1)
call phys_to_spectral_fg(a5f,gradphiy,1)
call phys_to_spectral_fg(a6f,gradphiz,1)

call dz_fg(gradphiz,a5)

!$acc parallel loop collapse(3)
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      a5(i,k,j,1)=a5(i,k,j,1)-kxpsi(i+cstartpsi(1))*gradphix(i,k,j,2) &
                           & -kypsi(j+cstartpsi(3))*gradphiy(i,k,j,2)
      a5(i,k,j,2)=a5(i,k,j,2)+kxpsi(i+cstartpsi(1))*gradphix(i,k,j,1) &
                           & +kypsi(j+cstartpsi(3))*gradphiy(i,k,j,1)
    enddo
  enddo
enddo

! assemble z component of surface tension force
!$acc parallel loop collapse(3)
do j=1,fpypsi
  do k=1,fpzpsi
    do i=1,npsix
      a4f(i,k,j)=-sigma(i,k,j)*fgradphix(i,k,j)*fgradphiz(i,k,j)
      a5f(i,k,j)=-sigma(i,k,j)*fgradphiy(i,k,j)*fgradphiz(i,k,j)
      a6f(i,k,j)=sigma(i,k,j)*((fgradphix(i,k,j))**2+(fgradphiy(i,k,j))**2)
    enddo
  enddo
enddo

call phys_to_spectral_fg(a4f,gradphix,1)
call phys_to_spectral_fg(a5f,gradphiy,1)
call phys_to_spectral_fg(a6f,gradphiz,1)

deallocate(fgradphix)
deallocate(fgradphiy)
deallocate(fgradphiz)
deallocate(sigma)
deallocate(a4f)
deallocate(a5f)
deallocate(a6f)

call dz_fg(gradphiz,a6)

!$acc parallel loop collapse(3)
do j=1,spypsi
  do k=1,npsiz
    do i=1,spxpsi
      a6(i,k,j,1)=a6(i,k,j,1)-kxpsi(i+cstartpsi(1))*gradphix(i,k,j,2) &
                           & -kypsi(j+cstartpsi(3))*gradphiy(i,k,j,2)
      a6(i,k,j,2)=a6(i,k,j,2)+kxpsi(i+cstartpsi(1))*gradphix(i,k,j,1) &
                           & +kypsi(j+cstartpsi(3))*gradphiy(i,k,j,1)
    enddo
  enddo
enddo


deallocate(gradphix)
deallocate(gradphiy)
deallocate(gradphiz)
! a4, a5, a6 on fine grid
! allocate variable in coarse grid for coupling with NS
allocate(gradphix(spx,nz,spy,2))
allocate(gradphiy(spx,nz,spy,2))
allocate(gradphiz(spx,nz,spy,2))

! take a4,a5,a6 from fine to coarse grid, then add to NS S term
call fine2coarse(a4,gradphix)
call fine2coarse(a5,gradphiy)
call fine2coarse(a6,gradphiz)

deallocate(a4)
deallocate(a5)
deallocate(a6)

#if match_dens == 2
! rescale NS equation if rhor > 1 for improved stability
!$acc kernels
s1=s1+3.0d0/dsqrt(8.0d0)*Ch/(we*rhor)*gradphix
s2=s2+3.0d0/dsqrt(8.0d0)*Ch/(we*rhor)*gradphiy
s3=s3+3.0d0/dsqrt(8.0d0)*Ch/(we*rhor)*gradphiz
!$acc end kernels
#else
!$acc kernels
s1=s1+3.0d0/dsqrt(8.0d0)*Ch/(we)*gradphix
s2=s2+3.0d0/dsqrt(8.0d0)*Ch/(we)*gradphiy
s3=s3+3.0d0/dsqrt(8.0d0)*Ch/(we)*gradphiz
!$acc end kernels
#endif

deallocate(gradphix)
deallocate(gradphiy)
deallocate(gradphiz)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! electric force calculation
#if eleflag == 1

allocate(gradphix(spx,nz,spy,2))
allocate(gradphiy(spx,nz,spy,2))
allocate(gradphiz(spx,nz,spy,2))
allocate(fgradphix(nx,fpz,fpy))
allocate(fgradphiy(nx,fpz,fpy))
allocate(fgradphiz(nx,fpz,fpy))

!$acc kernels
do i=1,spx
  gradphix(i,:,:,1)=-kx(i+cstart(1))*phic(i,:,:,2)
  gradphix(i,:,:,2)=+kx(i+cstart(1))*phic(i,:,:,1)
enddo
do j=1,spy
  gradphiy(:,:,j,1)=-ky(j+cstart(3))*phic(:,:,j,2)
  gradphiy(:,:,j,2)=+ky(j+cstart(3))*phic(:,:,j,1)
enddo
!$acc end kernels
call dz(phic,gradphiz)

call spectral_to_phys(gradphix,fgradphix,0)
call spectral_to_phys(gradphiy,fgradphiy,0)
call spectral_to_phys(gradphiz,fgradphiz,0)

!$acc kernels
threshold=0.9d0
! gthreshold=0.5d0*(1.0d0-threshold**2)/(dsqrt(2.0d0)*ch)
do j=1,fpy
 do k=1,fpz
  do i=1,nx
   modg=dabs(dsqrt((fgradphix(i,k,j))**2+(fgradphiy(i,k,j))**2+(fgradphiz(i,k,j))**2))
   gthreshold=0.75d0*(1.0d0-(phi(i,k,j))**2)/(dsqrt(2.0d0)*ch)
   if((dabs(phi(i,k,j)).lt.threshold).and.(modg.lt.gthreshold))then
    fgradphix(i,k,j)=fgradphix(i,k,j)/modg
    fgradphiy(i,k,j)=fgradphiy(i,k,j)/modg
    fgradphiz(i,k,j)=fgradphiz(i,k,j)/modg
   else
    fgradphix(i,k,j)=0.0d0
    fgradphiy(i,k,j)=0.0d0
    fgradphiz(i,k,j)=0.0d0
   endif
  enddo
 enddo
enddo
!$acc end kernels

call phys_to_spectral(fgradphix,gradphix,0)
call phys_to_spectral(fgradphiy,gradphiy,0)
call phys_to_spectral(fgradphiz,gradphiz,0)

deallocate(fgradphix,fgradphiy,fgradphiz)

#if match_dens == 2
! rescale NS equation if rhor > 1 for improved stability
!$acc kernels
s1=s1+stuart/rhor*gradphix
s2=s2+stuart/rhor*gradphiy
s3=s3+stuart/rhor*gradphiz
!$acc end kernels
#else
!$acc kernels
s1=s1+stuart*gradphix
s2=s2+stuart*gradphiy
s3=s3+stuart*gradphiz
!$acc end kernels
#endif

deallocate(gradphix,gradphiy,gradphiz)



! allocate(a4f(nx,fpz,fpy))
! allocate(a4(spx,nz,spy,2))
! ! charge density distribution, in [0,1]
! a4f=max(min(phi+1.0d0,1.0d0),0.0d0)
!
! call phys_to_spectral(a4f,a4,0)
!
! allocate(gradphix(spx,nz,spy,2))
! allocate(gradphiy(spx,nz,spy,2))
! allocate(gradphiz(spx,nz,spy,2))
! allocate(fgradphix(nx,fpz,fpy))
! allocate(fgradphiy(nx,fpz,fpy))
! allocate(fgradphiz(nx,fpz,fpy))
!
! do i=1,spx
!   gradphix(i,:,:,1)=-kx(i+cstart(1))*a4(i,:,:,2)
!   gradphix(i,:,:,2)=+kx(i+cstart(1))*a4(i,:,:,1)
! enddo
! do j=1,spy
!   gradphiy(:,:,j,1)=-ky(j+cstart(3))*a4(:,:,j,2)
!   gradphiy(:,:,j,2)=+ky(j+cstart(3))*a4(:,:,j,1)
! enddo
! call dz(a4,gradphiz)
!
! call spectral_to_phys(gradphix,fgradphix,0)
! call spectral_to_phys(gradphiy,fgradphiy,0)
! call spectral_to_phys(gradphiz,fgradphiz,0)
!
! fgradphix=a4f*fgradphix
! fgradphiy=a4f*fgradphiy
! fgradphiz=a4f*fgradphiz
!
! call phys_to_spectral(fgradphix,gradphix,0)
! call phys_to_spectral(fgradphiy,gradphiy,0)
! call phys_to_spectral(fgradphiz,gradphiz,0)
!
! deallocate(fgradphix,fgradphiy,fgradphiz)
! deallocate(a4,a4f)
!
! #if match_dens == 2
! ! rescale NS equation if rhor > 1 for improved stability
! s1=s1+stuart/rhor*gradphix
! s2=s2+stuart/rhor*gradphiy
! s3=s3+stuart/rhor*gradphiz
! #else
! s1=s1+stuart*gradphix
! s2=s2+stuart*gradphiy
! s3=s3+stuart*gradphiz
! #endif
!
! deallocate(gradphix,gradphiy,gradphiz)

#endif
! end of electric force calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! only for non-matched viscosities
! calculate non-linear part of viscous term (only for non-matched viscosities)
#if match_visc == 0
#if non_new_flag == 0
! if visr < 1
if (rank.eq.0) print*,'Unmatched viscosities - N/N'
! first row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call dz(uc,a6)

!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*uc(i,:,j,2)
    a4(i,:,j,2)=+kx(i+cstart(1))*uc(i,:,j,1)
    a5(i,:,j,1)=-ky(j+cstart(3))*uc(i,:,j,2)-kx(i+cstart(1))*vc(i,:,j,2)
    a5(i,:,j,2)=+ky(j+cstart(3))*uc(i,:,j,1)+kx(i+cstart(1))*vc(i,:,j,1)
    a6(i,:,j,1)=a6(i,:,j,1)-kx(i+cstart(1))*wc(i,:,j,2)
    a6(i,:,j,2)=a6(i,:,j,2)+kx(i+cstart(1))*wc(i,:,j,1)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=2.0d0*0.5d0*(visr-1.0d0)*(phif+1.0d0)*a4f(i,k,j)
      a5f(i,k,j)=0.5d0*(visr-1.0d0)*(phif+1.0d0)*a5f(i,k,j)
      a6f(i,k,j)=0.5d0*(visr-1.0d0)*(phif+1.0d0)*a6f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels

call phys_to_spectral(a4f,a4,0)
call phys_to_spectral(a5f,a5,0)
call phys_to_spectral(a6f,a6,0)

allocate(a7(spx,nz,spy,2))

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels

call dz(a6,a4)

! add to S term a4+a7
#if match_dens == 2
! rescale NS if rhor > 1 for improved stability
s1=s1+(a4+a7)/(re*rhor)
#else
s1=s1+(a4+a7)/re
#endif

deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(a7)

deallocate(a4f)
deallocate(a5f)
deallocate(a6f)


! second row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call dz(vc,a6)

!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-ky(j+cstart(3))*uc(i,:,j,2)-kx(i+cstart(1))*vc(i,:,j,2)
    a4(i,:,j,2)=+ky(j+cstart(3))*uc(i,:,j,1)+kx(i+cstart(1))*vc(i,:,j,1)
    a5(i,:,j,1)=-ky(j+cstart(3))*vc(i,:,j,2)
    a5(i,:,j,2)=+ky(j+cstart(3))*vc(i,:,j,1)
    a6(i,:,j,1)=a6(i,:,j,1)-ky(j+cstart(3))*wc(i,:,j,2)
    a6(i,:,j,2)=a6(i,:,j,2)+ky(j+cstart(3))*wc(i,:,j,1)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=0.5d0*(visr-1.0d0)*(phif+1.0d0)*a4f(i,k,j)
      a5f(i,k,j)=2.0d0*0.5d0*(visr-1.0d0)*(phif+1.0d0)*a5f(i,k,j)
      a6f(i,k,j)=0.5d0*(visr-1.0d0)*(phif+1.0d0)*a6f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels


call phys_to_spectral(a4f,a4,0)
call phys_to_spectral(a5f,a5,0)
call phys_to_spectral(a6f,a6,0)

allocate(a7(spx,nz,spy,2))

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels


call dz(a6,a4)

! add to S term a4+a7
#if match_dens == 2
s2=s2+(a4+a7)/(re*rhor)
#else
s2=s2+(a4+a7)/re
#endif

deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(a7)

deallocate(a4f)
deallocate(a5f)
deallocate(a6f)


! third row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call dz(uc,a4)
call dz(vc,a5)
call dz(wc,a6)

!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=a4(i,:,j,1)-kx(i+cstart(1))*wc(i,:,j,2)
    a4(i,:,j,2)=a4(i,:,j,2)+kx(i+cstart(1))*wc(i,:,j,1)
    a5(i,:,j,1)=a5(i,:,j,1)-ky(j+cstart(3))*wc(i,:,j,2)
    a5(i,:,j,2)=a5(i,:,j,2)+ky(j+cstart(3))*wc(i,:,j,1)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=0.5d0*(visr-1.0d0)*(phif+1.0d0)*a4f(i,k,j)
      a5f(i,k,j)=0.5d0*(visr-1.0d0)*(phif+1.0d0)*a5f(i,k,j)
      a6f(i,k,j)=2.0d0*0.5d0*(visr-1.0d0)*(phif+1.0d0)*a6f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels


call phys_to_spectral(a4f,a4,0)
call phys_to_spectral(a5f,a5,0)
call phys_to_spectral(a6f,a6,0)

allocate(a7(spx,nz,spy,2))

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels


call dz(a6,a4)

! add to S term a4+a7
#if match_dens == 2
!$acc kernels
s3=s3+(a4+a7)/(re*rhor)
!$acc end kernels
#else
!$acc kernels
s3=s3+(a4+a7)/re
!$acc end kernels
#endif

deallocate(a4,a5,a6,a7)
deallocate(a4f,a5f,a6f)

!End Block Visc-Contarst Newtonian
#endif


#if non_new_flag == 1
if (rank.eq.0) print*,'Unmatched viscosities - N/NN'

! first row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))
allocate(a7(spx,nz,spy,2))
allocate(a8(spx,nz,spy,2))
allocate(a9(spx,nz,spy,2))

call dz(uc,a7)
call dz(vc,a8)
call dz(wc,a9)


!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*2d0*uc(i,:,j,2)
    a4(i,:,j,2)=+kx(i+cstart(1))*2d0*uc(i,:,j,1)
    a5(i,:,j,1)=-ky(j+cstart(3))*uc(i,:,j,2)-kx(i+cstart(1))*vc(i,:,j,2)
    a5(i,:,j,2)=+ky(j+cstart(3))*uc(i,:,j,1)+kx(i+cstart(1))*vc(i,:,j,1)
    a6(i,:,j,1)=-ky(j+cstart(3))*2d0*vc(i,:,j,2)
    a6(i,:,j,2)=+ky(j+cstart(3))*2d0*vc(i,:,j,1)
    a7(i,:,j,1)=a7(i,:,j,1)-kx(i+cstart(1))*wc(i,:,j,2)
    a7(i,:,j,2)=a7(i,:,j,2)+kx(i+cstart(1))*wc(i,:,j,1)
    a8(i,:,j,1)=a8(i,:,j,1)-ky(j+cstart(3))*wc(i,:,j,2)
    a8(i,:,j,2)=a8(i,:,j,2)+ky(j+cstart(3))*wc(i,:,j,1)
    a9(i,:,j,1)=2d0*a9(i,:,j,1)
    a9(i,:,j,2)=2d0*a9(i,:,j,2)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))
allocate(a7f(nx,fpz,fpy))
allocate(a8f(nx,fpz,fpy))
allocate(a9f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)
call spectral_to_phys(a7,a7f,0)
call spectral_to_phys(a8,a8f,0)
call spectral_to_phys(a9,a9f,0)

deallocate(a4,a5,a6,a7,a8,a9)
allocate(viscnon(nx,fpz,fpy))

!computing the viscosity map
!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1d0,max(-1d0,phi(i,k,j)))
      gammadot=a4f(i,k,j)**2d0+a6f(i,k,j)**2d0+a9f(i,k,j)**2d0+&
      &          2d0*a5f(i,k,j)**2d0+2d0*a7f(i,k,j)**2d0+2d0*a8f(i,k,j)**2d0
      ! computing the equivalent shear rate...
      gammadot=dsqrt(0.5d0*gammadot)
      ! computing the dimensionelles viscosity
      viscnon(i,k,j)=muinfmuzero+(1d0-muinfmuzero)*((1d0+gammadot**2d0)**(0.5d0*(exp_non_new-1d0)))
    enddo
  enddo
enddo
!$acc end kernels

! computing the first row...
allocate(a10f(nx,fpz,fpy))
allocate(a11f(nx,fpz,fpy))
allocate(a12f(nx,fpz,fpy))

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1d0,max(-1d0,phi(i,k,j)))
      a10f(i,k,j)=0.5d0*(1d0+phif)*(viscnon(i,k,j)-1d0)*a4f(i,k,j)
      a11f(i,k,j)=0.5d0*(1d0+phif)*(viscnon(i,k,j)-1d0)*a5f(i,k,j)
      a12f(i,k,j)=0.5d0*(1d0+phif)*(viscnon(i,k,j)-1d0)*a7f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels


deallocate(a4f)

allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))
allocate(a7(spx,nz,spy,2))

call phys_to_spectral(a10f,a4,0)
call phys_to_spectral(a11f,a5,0)
call phys_to_spectral(a12f,a6,0)

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels


call dz(a6,a4)

! add to S term a4+a7
! rescale NS if rhor > 1 for improved stability
#if match_dens == 2
!$acc kernels
s1=s1+1.0d0/(re*rhor)*(a4+a7)
!$acc end kernels
#else
!$acc kernels
s1=s1+1.0d0/(re)*(a4+a7)
!$acc end kernels
#endif


! second row
!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1d0,max(-1d0,phi(i,k,j)))
      a10f(i,k,j)=0.5d0*(1.0d0+phif)*(viscnon(i,k,j)-1d0)*a5f(i,k,j)
      a11f(i,k,j)=0.5d0*(1.0d0+phif)*(viscnon(i,k,j)-1d0)*a6f(i,k,j)
      a12f(i,k,j)=0.5d0*(1.0d0+phif)*(viscnon(i,k,j)-1d0)*a8f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels


deallocate(a5f,a6f)

call phys_to_spectral(a10f,a4,0)
call phys_to_spectral(a11f,a5,0)
call phys_to_spectral(a12f,a6,0)

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels

call dz(a6,a4)

! add to S term a4+a7
! rescale NS if rhor > 1 for improved stability
#if match_dens == 2
!$acc kernels
s2=s2+1d0/(re*rhor)*(a4+a7)
!$acc end kernels
#else
!$acc kernels
s2=s2+1d0/(re)*(a4+a7)
!$acc end kernels
#endif

! third row
!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1d0,max(-1d0,phi(i,k,j)))
      a10f(i,k,j)=0.5d0*(1d0+phif)*(viscnon(i,k,j)-1d0)*a7f(i,k,j)
      a11f(i,k,j)=0.5d0*(1d0+phif)*(viscnon(i,k,j)-1d0)*a8f(i,k,j)
      a12f(i,k,j)=0.5d0*(1d0+phif)*(viscnon(i,k,j)-1d0)*a9f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels

deallocate(viscnon,a7f,a8f,a9f)

call phys_to_spectral(a10f,a4,0)
call phys_to_spectral(a11f,a5,0)
call phys_to_spectral(a12f,a6,0)

deallocate(a10f,a11f,a12f)

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels

call dz(a6,a4)

! add to S term a4+a7
! rescale NS if rhor > 1 for improved stability
#if match_dens == 2
!$acc kernels
s3=s3+1d0/(re*rhor)*(a4+a7)
!$acc end kernels
#else
!$acc kernels
s3=s3+1d0/(re)*(a4+a7)
!$acc end kernels
#endif


deallocate(a4,a5,a6,a7)
!Closing block non-newtonian
#endif







#elif match_visc == 2
! if visr > 1

! first row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call dz(uc,a6)

!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-kx(i+cstart(1))*uc(i,:,j,2)
    a4(i,:,j,2)=+kx(i+cstart(1))*uc(i,:,j,1)
    a5(i,:,j,1)=-ky(j+cstart(3))*uc(i,:,j,2)-kx(i+cstart(1))*vc(i,:,j,2)
    a5(i,:,j,2)=+ky(j+cstart(3))*uc(i,:,j,1)+kx(i+cstart(1))*vc(i,:,j,1)
    a6(i,:,j,1)=a6(i,:,j,1)-kx(i+cstart(1))*wc(i,:,j,2)
    a6(i,:,j,2)=a6(i,:,j,2)+kx(i+cstart(1))*wc(i,:,j,1)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=2.0d0*0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a4f(i,k,j)
      a5f(i,k,j)=0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a5f(i,k,j)
      a6f(i,k,j)=0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a6f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels

call phys_to_spectral(a4f,a4,0)
call phys_to_spectral(a5f,a5,0)
call phys_to_spectral(a6f,a6,0)

allocate(a7(spx,nz,spy,2))

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels

call dz(a6,a4)

! add to S term a4+a7
#if match_dens == 2
! rescale NS if rhor > 1 for improved stability
!$acc kernels
s1=s1+(a4+a7)/(re*rhor)*visr
!$acc end kernels
#else
!$acc kernels
s1=s1+(a4+a7)/re*visr
!$acc end kernels
#endif

deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(a7)

deallocate(a4f)
deallocate(a5f)
deallocate(a6f)


! second row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call dz(vc,a6)

!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=-ky(j+cstart(3))*uc(i,:,j,2)-kx(i+cstart(1))*vc(i,:,j,2)
    a4(i,:,j,2)=+ky(j+cstart(3))*uc(i,:,j,1)+kx(i+cstart(1))*vc(i,:,j,1)
    a5(i,:,j,1)=-ky(j+cstart(3))*vc(i,:,j,2)
    a5(i,:,j,2)=+ky(j+cstart(3))*vc(i,:,j,1)
    a6(i,:,j,1)=a6(i,:,j,1)-ky(j+cstart(3))*wc(i,:,j,2)
    a6(i,:,j,2)=a6(i,:,j,2)+ky(j+cstart(3))*wc(i,:,j,1)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a4f(i,k,j)
      a5f(i,k,j)=2.0d0*0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a5f(i,k,j)
      a6f(i,k,j)=0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a6f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels

call phys_to_spectral(a4f,a4,0)
call phys_to_spectral(a5f,a5,0)
call phys_to_spectral(a6f,a6,0)

allocate(a7(spx,nz,spy,2))

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels

call dz(a6,a4)

! add to S term a4+a7
#if match_dens == 2
! rescale NS if rhor > 1 for improved stability
!$acc kernels
s2=s2+(a4+a7)/(re*rhor)*visr
!$acc end kernels
#else
!$acc kernels
s2=s2+(a4+a7)/re*visr
!$acc end kernels
#endif

deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(a7)

deallocate(a4f)
deallocate(a5f)
deallocate(a6f)


! third row
allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call dz(uc,a4)
call dz(vc,a5)
call dz(wc,a6)

!$acc kernels
do j=1,spy
  do i=1,spx
    a4(i,:,j,1)=a4(i,:,j,1)-kx(i+cstart(1))*wc(i,:,j,2)
    a4(i,:,j,2)=a4(i,:,j,2)+kx(i+cstart(1))*wc(i,:,j,1)
    a5(i,:,j,1)=a5(i,:,j,1)-ky(j+cstart(3))*wc(i,:,j,2)
    a5(i,:,j,2)=a5(i,:,j,2)+ky(j+cstart(3))*wc(i,:,j,1)
  enddo
enddo
!$acc end kernels


allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(a4,a4f,0)
call spectral_to_phys(a5,a5f,0)
call spectral_to_phys(a6,a6f,0)

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a4f(i,k,j)
      a5f(i,k,j)=0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a5f(i,k,j)
      a6f(i,k,j)=2.0d0*0.5d0*(1.0d0-1.0d0/visr)*(phif-1.0d0)*a6f(i,k,j)
    enddo
  enddo
enddo
!$acc end kernels


call phys_to_spectral(a4f,a4,0)
call phys_to_spectral(a5f,a5,0)
call phys_to_spectral(a6f,a6,0)

allocate(a7(spx,nz,spy,2))

!$acc kernels
do j=1,spy
  do i=1,spx
    a7(i,:,j,1)=-kx(i+cstart(1))*a4(i,:,j,2)-ky(j+cstart(3))*a5(i,:,j,2)
    a7(i,:,j,2)=+kx(i+cstart(1))*a4(i,:,j,1)+ky(j+cstart(3))*a5(i,:,j,1)
  enddo
enddo
!$acc end kernels


call dz(a6,a4)

! add to S term a4+a7
#if match_dens == 2
! rescale NS if rhor > 1 for improved stability
!$acc kernels
s3=s3+(a4+a7)/(re*rhor)*visr
!$acc end kernels
#else
!$acc kernels
s3=s3+(a4+a7)/re*visr
!$acc end kernels
#endif

deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(a7)

deallocate(a4f)
deallocate(a5f)
deallocate(a6f)

#endif


! only for non-matched densities
#if match_dens != 1
! calculate gravity and buoyancy terms (only for non-matched densities)
#if match_dens == 0
#if b_type == 0
! no gravity, no further terms added to S

#elif b_type == 1
! gravity and buoyancy
! gravity array is [x,z,y] to keep the array ordering as usual in the code
! S term order is S1:x, S2:y, S3:z
allocate(a4f(nx,fpz,fpy))

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=1.0d0+0.5d0*(rhor-1.0d0)*(phif+1.0d0)
    enddo
  enddo
enddo
!$acc end kernels


allocate(a4(spx,nz,spy,2))

call phys_to_spectral(a4f,a4,0)

deallocate(a4f)

!$acc kernels
s1=s1+1.0d0/fr**2*a4*grav(1)
s2=s2+1.0d0/fr**2*a4*grav(3)
s3=s3+1.0d0/fr**2*a4*grav(2)
!$acc end kernels


deallocate(a4)

#elif b_type == 2
! only buoyancy
! gravity array is [x,z,y] to keep the array ordering as usual in the code
! S term order is S1:x, S2:y, S3:z
allocate(a4f(nx,fpz,fpy))

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=0.5d0*(rhor-1.0d0)*(phif+1.0d0)
    enddo
  enddo
enddo
!$acc end kernels


allocate(a4(spx,nz,spy,2))

call phys_to_spectral(a4f,a4,0)

deallocate(a4f)

!$acc kernels
s1=s1+1.0d0/fr**2*a4*grav(1)
s2=s2+1.0d0/fr**2*a4*grav(3)
s3=s3+1.0d0/fr**2*a4*grav(2)
!$acc end kernels


deallocate(a4)

#endif
#elif match_dens == 2
! rescale NS if rhor > 1 for improved stability
#if b_type == 0
! no gravity, no further terms added to S

#elif b_type == 1
! gravity and buoyancy
! gravity array is [x,z,y] to keep the array ordering as usual in the code
! S term order is S1:x, S2:y, S3:z
allocate(a4f(nx,fpz,fpy))

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=1.0d0+0.5d0*(1.0d0/rhor-1.0d0)*(1.0d0-phif)
    enddo
  enddo
enddo
!$acc end kernels


allocate(a4(spx,nz,spy,2))

call phys_to_spectral(a4f,a4,0)

deallocate(a4f)

!$acc kernels
s1=s1+1.0d0/fr**2*a4*grav(1)
s2=s2+1.0d0/fr**2*a4*grav(3)
s3=s3+1.0d0/fr**2*a4*grav(2)
!$acc end kernels


deallocate(a4)

#elif b_type == 2
! only buoyancy
! gravity array is [x,z,y] to keep the array ordering as usual in the code
! S term order is S1:x, S2:y, S3:z
allocate(a4f(nx,fpz,fpy))

!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
! delta rho with respect to carrier: rho-rho_c
      a4f(i,k,j)=0.5d0*(1.0d0-1.0d0/rhor)*(1.0d0+phif)
    enddo
  enddo
enddo
!$acc end kernels


allocate(a4(spx,nz,spy,2))

call phys_to_spectral(a4f,a4,0)

deallocate(a4f)

s1=s1+1.0d0/fr**2*a4*grav(1)
s2=s2+1.0d0/fr**2*a4*grav(3)
s3=s3+1.0d0/fr**2*a4*grav(2)

deallocate(a4)

#endif
#endif

! calculate non-linear part of the time derivative term (only for non-matched densities)

! phi already available in physical space from surface force calculation, no need to transform it again
allocate(a4f(nx,fpz,fpy))
allocate(a5f(nx,fpz,fpy))
allocate(a6f(nx,fpz,fpy))

call spectral_to_phys(uc-ucp,a4f,1)
call spectral_to_phys(vc-vcp,a5f,1)
call spectral_to_phys(wc-wcp,a6f,1)

#if match_dens == 0
!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=(phif+1.0d0)*(a4f(i,k,j))/dt
      a5f(i,k,j)=(phif+1.0d0)*(a5f(i,k,j))/dt
      a6f(i,k,j)=(phif+1.0d0)*(a6f(i,k,j))/dt
    enddo
  enddo
enddo
!$acc end kernels
#elif match_dens == 2
!$acc kernels
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=(1.0d0-phif)*(a4f(i,k,j))/dt
      a5f(i,k,j)=(1.0d0-phif)*(a5f(i,k,j))/dt
      a6f(i,k,j)=(1.0d0-phif)*(a6f(i,k,j))/dt
    enddo
  enddo
enddo
!$acc end kernels
#endif


allocate(a4(spx,nz,spy,2))
allocate(a5(spx,nz,spy,2))
allocate(a6(spx,nz,spy,2))

call phys_to_spectral(a4f,a4,1)
call phys_to_spectral(a5f,a5,1)
call phys_to_spectral(a6f,a6,1)

deallocate(a4f)
deallocate(a5f)
deallocate(a6f)

#if match_dens == 2
! rescale NS if rhor > 1 for improved stability
!$acc kernels
s1=s1-0.5d0*(1.0d0/rhor-1.0d0)*a4
s2=s2-0.5d0*(1.0d0/rhor-1.0d0)*a5
s3=s3-0.5d0*(1.0d0/rhor-1.0d0)*a6
!$acc end kernels
#else
!$acc kernels
s1=s1-0.5d0*(rhor-1.0d0)*a4
s2=s2-0.5d0*(rhor-1.0d0)*a5
s3=s3-0.5d0*(rhor-1.0d0)*a6
!$acc end kernels
#endif

deallocate(a4)
deallocate(a5)
deallocate(a6)

! update velocity previous value; uc,vc,wc will be updated in the following calculation
!$acc kernels
ucp=uc
vcp=vc
wcp=wc
!$acc end kernels

#endif


! add body force contribution to N-S non-linear terms (kind of Boussinesq approximation)
! body force direction array is [x,z,y] to keep the array ordering as usual in the code
! only active in phi=+1
#if bodyflag == 1
allocate(a4(spx,nz,spy,2))
call phys_to_spectral(phi+1.0d0,a4,0)
!$acc kernels
s1=s1+body_d(1)*body_c*0.5d0*a4
s2=s2+body_d(3)*body_c*0.5d0*a4
s3=s3+body_d(2)*body_c*0.5d0*a4
!$acc end kernels
deallocate(a4)
#endif


! add s-shaped pressure gradient contribution to N-S non-linear terms
! body force direction array is [x,z,y] to keep the array ordering as usual in the code
! no sense to have pressure gradient along z, i.e. sgradp_d
#if sgradpflag == 1
allocate(a4(spx,nz,spy,2))
call phys_to_spectral(phi,a4,0)
!$acc kernels
s1=s1+sgradp_d(1)*a4
s2=s2+sgradp_d(3)*a4
!$acc end kernels
deallocate(a4)
#endif

return
end
