subroutine get_sterm

use commondata
use velocity
use wavenumbers
use sterm

double precision, dimension(nx/2+1,nz,ny,2) :: s1s,s2s,s3s

integer :: i,j,k

! pressure gradient
do j=1,ny
  do k=1,nz
    do i=1,nx/2+1
      s1s(i,k,j,1)=-kx(i)*pressc(i,k,j,2)
      s1s(i,k,j,2)=kx(i)*pressc(i,k,j,1)
      s2s(i,k,j,1)=-ky(j)*pressc(i,k,j,2)
      s2s(i,k,j,2)=ky(j)*pressc(i,k,j,1)
    enddo
  enddo
enddo

call dz(pressc,s3s)

s1s(1,1,1,1)=-s1s(1,1,1,1)-gradpx*dble(2*nx*ny)
s2s(1,1,1,1)=-s2s(1,1,1,1)-gradpy*dble(2*nx*ny)
s3s=-s3s

! viscous term, surface force, gravity
call phi_non_linear(s1s,s2s,s3s)


call spectral_to_phys(s1s,s1,0)
call spectral_to_phys(s2s,s2,0)
call spectral_to_phys(s3s,s3,0)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phi_non_linear(s1,s2,s3)

use commondata
use wavenumbers
use velocity

double precision, dimension(nx/2+1,nz,ny,2) :: s1,s2,s3,one_s
double precision, allocatable, dimension(:,:,:,:) :: gradphix,gradphiy,gradphiz,nablaphi
double precision, allocatable, dimension(:,:,:,:) :: a4,a5,a6,a7,a8
double precision, allocatable, dimension(:,:,:) :: fgradphix,fgradphiy,fgradphiz,fnablaphi
double precision, allocatable, dimension(:,:,:) :: a4f,a5f,a6f,a7f,a8f
double precision, allocatable, dimension(:,:,:,:) :: gradfunx,gradfuny,gradfunz,func
double precision, allocatable, dimension(:,:,:) :: fgradfunx,fgradfuny,fgradfunz
double precision, allocatable, dimension(:,:,:)   :: fun,forcex,forcey,forcez


integer :: i,j,k,fpy,fpz,spx,spy,cstart(3),fstart(3)

fpy=ny
fpz=nz
spx=nx/2+1
spy=ny
cstart=0
fstart=0
one_s=0.0d0
one_s(1,1,1,1)=dble(2*nx*ny)

! calculate surface force
call spectral_to_phys(phic,phi,1)

allocate(nablaphi(spx,nz,spy,2))
allocate(gradphix(spx,nz,spy,2))


! use gradphix to store temporarily z derivatives of phi
call dz(phic,gradphix)
call dz(gradphix,nablaphi)
gradphix=0.0d0

do j=1,spy
  do i=1,spx
    nablaphi(i,:,j,1)=nablaphi(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,1)
    nablaphi(i,:,j,2)=nablaphi(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*phic(i,:,j,2)
  enddo
enddo


allocate(fnablaphi(nx,fpz,fpy))
call spectral_to_phys(nablaphi,fnablaphi,1)
deallocate(nablaphi)


allocate(gradphiy(spx,nz,spy,2))
! calculate the gradient of phi (x,y components)
do j=1,spy
  do i=1,spx
    gradphix(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
    gradphix(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
    gradphiy(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
    gradphiy(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
  enddo
enddo

allocate(fgradphix(nx,fpz,fpy))
call spectral_to_phys(gradphix,fgradphix,1)
deallocate(gradphix)

allocate(fgradphiy(nx,fpz,fpy))
call spectral_to_phys(gradphiy,fgradphiy,1)
deallocate(gradphiy)


! calculate the gradient of phi (z component)
allocate(gradphiz(spx,nz,spy,2))
call dz(phic,gradphiz)
allocate(fgradphiz(nx,fpz,fpy))
call spectral_to_phys(gradphiz,fgradphiz,1)
deallocate(gradphiz)


if(psiflag.eq.0)then
  ! calculate chemical potential times the gradient of phi and update the gradient of phi
  ! (phi**3-phi-Ch**2*nablaphi)*grad phi
  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        fgradphix(i,k,j)=(phi(i,k,j)**3-phi(i,k,j)-ch**2*fnablaphi(i,k,j))*fgradphix(i,k,j)
        fgradphiy(i,k,j)=(phi(i,k,j)**3-phi(i,k,j)-ch**2*fnablaphi(i,k,j))*fgradphiy(i,k,j)
        fgradphiz(i,k,j)=(phi(i,k,j)**3-phi(i,k,j)-ch**2*fnablaphi(i,k,j))*fgradphiz(i,k,j)
      enddo
    enddo
  enddo
  deallocate(fnablaphi)

  allocate(gradphix(spx,nz,spy,2))
  allocate(gradphiy(spx,nz,spy,2))
  allocate(gradphiz(spx,nz,spy,2))

  call phys_to_spectral(fgradphix,gradphix,1)
  call phys_to_spectral(fgradphiy,gradphiy,1)
  call phys_to_spectral(fgradphiz,gradphiz,1)
elseif(psiflag.eq.1)then
  call spectral_to_phys(psic,psi,1)
  allocate(fun(nx,fpz,fpy))
  fun=1.0d0-El*psi
  allocate(forcex(nx,fpz,fpy))
  allocate(forcey(nx,fpz,fpy))
  allocate(forcez(nx,fpz,fpy))
  !! CAPILLARY FORCE (1-\psi)/(We*Ch) NORMAL COMPONENT
  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        forcex(i,k,j)=fun(i,k,j)*(phi(i,k,j)**3-phi(i,k,j)-ch**2*fnablaphi(i,k,j))*fgradphix(i,k,j)
        forcey(i,k,j)=fun(i,k,j)*(phi(i,k,j)**3-phi(i,k,j)-ch**2*fnablaphi(i,k,j))*fgradphiy(i,k,j)
        forcez(i,k,j)=fun(i,k,j)*(phi(i,k,j)**3-phi(i,k,j)-ch**2*fnablaphi(i,k,j))*fgradphiz(i,k,j)
      enddo
    enddo
  enddo
  deallocate(fnablaphi)
  !! MARANGONI FORCE TANGENTIAL COMPONENT
  !! FUNCTION OF SURFACE TENSION (PSI)
  allocate(func(spx,nz,spy,2))

  call phys_to_spectral(fun,func,1)

  deallocate(fun)

  allocate(gradfunx(spx,nz,spy,2))
  allocate(gradfuny(spx,nz,spy,2))
  allocate(gradfunz(spx,nz,spy,2))
  !! COMPUTING THE GRADIENT OF fun.
  do j=1,spy
    do i=1,spx
      gradfunx(i,:,j,1)=-kx(i+cstart(1))*func(i,:,j,2)
      gradfunx(i,:,j,2)=kx(i+cstart(1))*func(i,:,j,1)
      gradfuny(i,:,j,1)=-ky(j+cstart(3))*func(i,:,j,2)
      gradfuny(i,:,j,2)=ky(j+cstart(3))*func(i,:,j,1)
    enddo
  enddo

  call dz(func,gradfunz)

  deallocate(func)

  allocate(fgradfunx(nx,fpz,fpy))
  allocate(fgradfuny(nx,fpz,fpy))
  allocate(fgradfunz(nx,fpz,fpy))

  call spectral_to_phys(gradfunx,fgradfunx,1)
  call spectral_to_phys(gradfuny,fgradfuny,1)
  call spectral_to_phys(gradfunz,fgradfunz,1)

  deallocate(gradfunx)
  deallocate(gradfuny)
  deallocate(gradfunz)

  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        forcex(i,k,j)=forcex(i,k,j)+ch**2*((fgradphiy(i,k,j)**2+fgradphiz(i,k,j)**2)*fgradfunx(i,k,j) &
  &      - (fgradphix(i,k,j)*fgradphiy(i,k,j))*fgradfuny(i,k,j) - (fgradphix(i,k,j)*fgradphiz(i,k,j))*fgradfunz(i,k,j))
        forcey(i,k,j)=forcey(i,k,j)+ch**2*((-fgradphix(i,k,j)*fgradphiy(i,k,j))*fgradfunx(i,k,j) &
  &      + (fgradphix(i,k,j)**2 + fgradphiz(i,k,j)**2)*fgradfuny(i,k,j) - (fgradphiy(i,k,j)*fgradphiz(i,k,j))*fgradfunz(i,k,j))
        forcez(i,k,j)=forcez(i,k,j)+ch**2*((-fgradphix(i,k,j)*fgradphiz(i,k,j))*fgradfunx(i,k,j) &
  &      - (fgradphiy(i,k,j)*fgradphiz(i,k,j))*fgradfuny(i,k,j) + (fgradphix(i,k,j)**2 + fgradphiy(i,k,j)**2)*fgradfunz(i,k,j))
      enddo
    enddo
  enddo

  allocate(gradphix(spx,nz,spy,2))
  allocate(gradphiy(spx,nz,spy,2))
  allocate(gradphiz(spx,nz,spy,2))

  call phys_to_spectral(forcex,gradphix,1)
  call phys_to_spectral(forcey,gradphiy,1)
  call phys_to_spectral(forcez,gradphiz,1)
endif





s1=s1+3.0d0/sqrt(8.0d0)*1.0d0/(we*ch)*gradphix
s2=s2+3.0d0/sqrt(8.0d0)*1.0d0/(we*ch)*gradphiy
s3=s3+3.0d0/sqrt(8.0d0)*1.0d0/(we*ch)*gradphiz

deallocate(gradphix)
deallocate(gradphiy)
deallocate(gradphiz)
deallocate(fgradphix)
deallocate(fgradphiy)
deallocate(fgradphiz)

! add viscous term:
allocate(gradphix(nx/2+1,nz,ny,2))
allocate(gradphiy(nx/2+1,nz,ny,2))
allocate(gradphiz(nx/2+1,nz,ny,2))
allocate(gradfunx(nx/2+1,nz,ny,2))
allocate(gradfuny(nx/2+1,nz,ny,2))
allocate(gradfunz(nx/2+1,nz,ny,2))

call dz(uc,gradphix)
call dz(gradphix,gradfunx)
deallocate(gradphix)
call dz(vc,gradphiy)
call dz(gradphiy,gradfuny)
deallocate(gradphiy)
call dz(wc,gradphiz)
call dz(gradphiz,gradfunz)
deallocate(gradphiz)

do j=1,ny
  do k=1,nz
    do i=1,nx/2+1
      s1(i,k,j,1)=s1(i,k,j,1)+(gradfunx(i,k,j,1)-k2(i,j)*uc(i,k,j,1))/Re
      s1(i,k,j,2)=s1(i,k,j,2)+(gradfunx(i,k,j,2)-k2(i,j)*uc(i,k,j,2))/Re
      s2(i,k,j,1)=s2(i,k,j,1)+(gradfuny(i,k,j,1)-k2(i,j)*vc(i,k,j,1))/Re
      s2(i,k,j,2)=s2(i,k,j,2)+(gradfuny(i,k,j,2)-k2(i,j)*vc(i,k,j,2))/Re
      s3(i,k,j,1)=s3(i,k,j,1)+(gradfunz(i,k,j,1)-k2(i,j)*wc(i,k,j,1))/Re
      s3(i,k,j,2)=s3(i,k,j,2)+(gradfunz(i,k,j,2)-k2(i,j)*wc(i,k,j,2))/Re
    enddo
  enddo
enddo

deallocate(gradfunx)
deallocate(gradfuny)
deallocate(gradfunz)

! only for non-matched viscosities
if(abs(visr-1.0d0).ge.1e-8)then
  ! calculate non-linear part of viscous term (only for non-matched viscosities)
  ! first part: (phi+1)*Nabla u
  allocate(gradphix(spx,nz,spy,2))
  allocate(gradphiy(spx,nz,spy,2))
  allocate(gradphiz(spx,nz,spy,2))
  allocate(nablaphi(spx,nz,spy,2))

  call dz(uc,nablaphi)
  call dz(nablaphi,gradphix)

  call dz(vc,nablaphi)
  call dz(nablaphi,gradphiy)

  call dz(wc,nablaphi)
  call dz(nablaphi,gradphiz)

  deallocate(nablaphi)

  do j=1,spy
    do i=1,spx
      gradphix(i,:,j,1)=gradphix(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*uc(i,:,j,1)
      gradphix(i,:,j,2)=gradphix(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*uc(i,:,j,2)
      gradphiy(i,:,j,1)=gradphiy(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*vc(i,:,j,1)
      gradphiy(i,:,j,2)=gradphiy(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*vc(i,:,j,2)
      gradphiz(i,:,j,1)=gradphiz(i,:,j,1)-k2(i+cstart(1),j+cstart(3))*wc(i,:,j,1)
      gradphiz(i,:,j,2)=gradphiz(i,:,j,2)-k2(i+cstart(1),j+cstart(3))*wc(i,:,j,2)
    enddo
  enddo

  allocate(fgradphix(nx,fpz,fpy))
  allocate(fgradphiy(nx,fpz,fpy))
  allocate(fgradphiz(nx,fpz,fpy))

  call spectral_to_phys(gradphix,fgradphix,1)
  call spectral_to_phys(gradphiy,fgradphiy,1)
  call spectral_to_phys(gradphiz,fgradphiz,1)

  ! phi already available in physical space from surface force calculation, no need to transform it again
  ! assemble first part of non-linear term
  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        fgradphix(i,k,j)=(phi(i,k,j)+1.0d0)*fgradphix(i,k,j)
        fgradphiy(i,k,j)=(phi(i,k,j)+1.0d0)*fgradphiy(i,k,j)
        fgradphiz(i,k,j)=(phi(i,k,j)+1.0d0)*fgradphiz(i,k,j)
      enddo
    enddo
  enddo

  call phys_to_spectral(fgradphix,gradphix,1)
  call phys_to_spectral(fgradphiy,gradphiy,1)
  call phys_to_spectral(fgradphiz,gradphiz,1)

  s1=s1+(visr-1.0d0)/(2.0d0*re)*gradphix
  s2=s2+(visr-1.0d0)/(2.0d0*re)*gradphiy
  s3=s3+(visr-1.0d0)/(2.0d0*re)*gradphiz

  deallocate(gradphix)
  deallocate(gradphiy)
  deallocate(gradphiz)

  deallocate(fgradphix)
  deallocate(fgradphiy)
  deallocate(fgradphiz)

  ! second part: grad phi * (grad u + grad u')
  allocate(gradphix(spx,nz,spy,2))
  allocate(gradphiy(spx,nz,spy,2))
  allocate(gradphiz(spx,nz,spy,2))

  ! calculate phi derivatives
  call dz(phic,gradphiz)
  do j=1,spy
    do i=1,spx
      gradphix(i,:,j,1)=-kx(i+cstart(1))*phic(i,:,j,2)
      gradphix(i,:,j,2)=kx(i+cstart(1))*phic(i,:,j,1)
      gradphiy(i,:,j,1)=-ky(j+cstart(3))*phic(i,:,j,2)
      gradphiy(i,:,j,2)=ky(j+cstart(3))*phic(i,:,j,1)
    enddo
  enddo

  allocate(fgradphix(nx,fpz,fpy))
  allocate(fgradphiy(nx,fpz,fpy))
  allocate(fgradphiz(nx,fpz,fpy))

  call spectral_to_phys(gradphix,fgradphix,1)
  call spectral_to_phys(gradphiy,fgradphiy,1)
  call spectral_to_phys(gradphiz,fgradphiz,1)

  deallocate(gradphix)
  deallocate(gradphiy)
  deallocate(gradphiz)


  allocate(a4(spx,nz,spy,2))
  allocate(a5(spx,nz,spy,2))
  allocate(a6(spx,nz,spy,2))
  allocate(a7(spx,nz,spy,2))
  allocate(a8(spx,nz,spy,2))

  call dz(uc,a6)
  call dz(vc,a8)
  do j=1,spy
    do i=1,spx
      a4(i,:,j,1)=-2.0d0*kx(i+cstart(1))*uc(i,:,j,2)
      a4(i,:,j,2)=2.0d0*kx(i+cstart(1))*uc(i,:,j,1)
      a5(i,:,j,1)=-ky(j+cstart(3))*uc(i,:,j,2)-kx(i+cstart(1))*vc(i,:,j,2)
      a5(i,:,j,2)=ky(j+cstart(3))*uc(i,:,j,1)+kx(i+cstart(1))*vc(i,:,j,1)
      a6(i,:,j,1)=a6(i,:,j,1)-kx(i+cstart(1))*wc(i,:,j,2)
      a6(i,:,j,2)=a6(i,:,j,2)+kx(i+cstart(1))*wc(i,:,j,1)
      a7(i,:,j,1)=-2.0d0*ky(j+cstart(3))*vc(i,:,j,2)
      a7(i,:,j,2)=2.0d0*ky(j+cstart(3))*vc(i,:,j,1)
      a8(i,:,j,1)=a8(i,:,j,1)-ky(j+cstart(3))*wc(i,:,j,2)
      a8(i,:,j,2)=a8(i,:,j,2)+ky(j+cstart(3))*wc(i,:,j,1)
    enddo
  enddo

  allocate(a4f(nx,fpz,fpy))
  allocate(a5f(nx,fpz,fpy))
  allocate(a6f(nx,fpz,fpy))
  allocate(a7f(nx,fpz,fpy))
  allocate(a8f(nx,fpz,fpy))

  call spectral_to_phys(a4,a4f,1)
  call spectral_to_phys(a7,a7f,1)
  call spectral_to_phys(a5,a5f,1)
  call spectral_to_phys(a6,a6f,1)
  call spectral_to_phys(a8,a8f,1)

  deallocate(a4)
  deallocate(a7)

  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        a4f(i,k,j)=fgradphix(i,k,j)*a4f(i,k,j)+fgradphiy(i,k,j)*a5f(i,k,j)+fgradphiz(i,k,j)*a6f(i,k,j)
        a7f(i,k,j)=fgradphix(i,k,j)*a5f(i,k,j)+fgradphiy(i,k,j)*a7f(i,k,j)+fgradphiz(i,k,j)*a8f(i,k,j)
      enddo
    enddo
  enddo


  call dz(wc,a5)
  a5=2.0d0*a5
  call spectral_to_phys(a5,a5f,1)

  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        a5f(i,k,j)=fgradphix(i,k,j)*a6f(i,k,j)+fgradphiy(i,k,j)*a8f(i,k,j)+fgradphiz(i,k,j)*a5f(i,k,j)
      enddo
    enddo
  enddo

  deallocate(a6f)
  deallocate(a8f)

  ! a4f: first row of matrix in physical space
  ! a7f: second row of matrix in physical space
  ! a5f: third row of matrix in physical space

  call phys_to_spectral(a4f,a5,1)
  call phys_to_spectral(a7f,a6,1)
  call phys_to_spectral(a5f,a8,1)

  ! a5: first row of matrix in spectral space
  ! a6: second row of matrix in spectral space
  ! a8: third row of matrix in spectral space

  deallocate(a4f)
  deallocate(a5f)
  deallocate(a7f)

  s1=s1+(visr-1.0d0)/(2.0d0*re)*a5
  s2=s2+(visr-1.0d0)/(2.0d0*re)*a6
  s3=s3+(visr-1.0d0)/(2.0d0*re)*a8

  deallocate(a5)
  deallocate(a6)
  deallocate(a8)
endif


! only for non-matched densities
if(abs(rhor-1.0d0).ge.1e-8)then
  ! calculate gravity and buoyancy terms (only for non-matched densities)
  if(b_type.eq.0)then
    ! no gravity, no further terms added to S
  elseif(b_type.eq.1)then
    ! gravity and buoyancy
    ! gravity array is [x,z,y] to keep the array ordering as usual in the code
    ! S term order is S1:x, S2:y, S3:z
    s1=s1+1.0d0/fr**2*(one_s+0.5d0*(rhor-1.0d0)*(phic+one_s))*grav(1)
    s2=s2+1.0d0/fr**2*(one_s+0.5d0*(rhor-1.0d0)*(phic+one_s))*grav(3)
    s3=s3+1.0d0/fr**2*(one_s+0.5d0*(rhor-1.0d0)*(phic+one_s))*grav(2)
  elseif(b_type.eq.2)then
    ! only buoyancy
    ! gravity array is [x,z,y] to keep the array ordering as usual in the code
    ! S term order is S1:x, S2:y, S3:z
    s1=s1+1.0d0/fr**2*(0.5d0*(rhor-1.0d0)*(phic+one_s))*grav(1)
    s2=s2+1.0d0/fr**2*(0.5d0*(rhor-1.0d0)*(phic+one_s))*grav(3)
    s3=s3+1.0d0/fr**2*(0.5d0*(rhor-1.0d0)*(phic+one_s))*grav(2)
  endif
endif




return
end
