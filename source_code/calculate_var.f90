subroutine calculate_w(h)

use commondata
use par_size
use wavenumber
use velocity
use sim_par
use velocity

double precision :: beta2(spx,spy),k2l(spx,spy)
double precision :: h(spx,nz,spy,2)
double precision, allocatable :: dw1(:,:,:,:),dw2(:,:,:),dw3(:,:,:)
double precision, allocatable :: ddw1(:,:,:,:),ddw2(:,:,:),ddw3(:,:,:)
double precision, allocatable, dimension(:,:) :: a11,a12,a21,a22,fr1,fr2,fc1,fc2
double precision, allocatable :: A(:,:,:),B(:,:,:)
double precision :: det

integer :: i,j,k

h=h/gamma

do j=1,spy
  do i=1,spx
    beta2(i,j)=(1.0d0+gamma*k2(i+cstart(1),j+cstart(3)))/gamma
    k2l(i,j)=k2(i+cstart(1),j+cstart(3))
  enddo
enddo

! solve Helmholtz equation for auxiliary variable
call helmholtz(h,beta2,p_w,q_w,r_w,zp)

! solve Helmholtz equation for w_1 (see notes)
call helmholtz(h,k2l,p_w,q_w,r_w,zp)


! calculate first derivative of w_1,w_2,w_3 along z
allocate(dw1(spx,nz,spy,2))
allocate(dw2(spx,nz,spy))
allocate(dw3(spx,nz,spy))
call dz(h,dw1)
call dz_red(wa2,dw2)
call dz_red(wa3,dw3)


! set boundary conditions, compile.sh takes care of substituting
! correct values for conditional compilation
#define bc_up boundary_conditions_up
#define bc_down boundary_conditions_down
! solve system:
!  _       _   _ _   _         _
! | a11 a12 |*| A |=| fr1+i*fc1 |
! |_a21 a22_| |_B_| |_fr2+i*fc1_|
!
! fr1,fr2 are the real part, fc1,fc2 the imaginary part, A,B are complex valued


allocate(a11(spx,spy))
allocate(a12(spx,spy))
allocate(a21(spx,spy))
allocate(a22(spx,spy))
allocate(fr1(spx,spy))
allocate(fr2(spx,spy))
allocate(fc1(spx,spy))
allocate(fc2(spx,spy))

allocate(ddw1(spx,nz,spy,2))
allocate(ddw2(spx,nz,spy))
allocate(ddw3(spx,nz,spy))

! lower boundary: z=-1
#if bc_down==0
! no-slip at z=-1
a11(:,:)=0.5d0*dw2(:,1,:)
a12(:,:)=0.5d0*dw3(:,1,:)
fr1(:,:)=0.5d0*dw1(:,1,:,1)
fc1(:,:)=0.5d0*dw1(:,1,:,2)
do k=2,nz
! -1**(k-1) needed to evaluate derivative at boundary, see Canuto et al. 2006, pag. 85
  a11(:,:)=a11(:,:)+(zp(1))**(k-1)*dw2(:,k,:)
  a12(:,:)=a12(:,:)+(zp(1))**(k-1)*dw3(:,k,:)
  fr1(:,:)=fr1(:,:)+(zp(1))**(k-1)*dw1(:,k,:,1)
  fc1(:,:)=fc1(:,:)+(zp(1))**(k-1)*dw1(:,k,:,2)
enddo
#else
! free-slip at z=-1
call dz(dw1,ddw1)
call dz_red(dw2,ddw2)
call dz_red(dw3,ddw3)
a11(:,:)=0.5d0*ddw2(:,1,:)
a12(:,:)=0.5d0*ddw3(:,1,:)
fr1(:,:)=0.5d0*ddw1(:,1,:,1)
fc1(:,:)=0.5d0*ddw1(:,1,:,2)
do k=2,nz
  a11(:,:)=a11(:,:)+(zp(1))**(k-1)*ddw2(:,k,:)
  a12(:,:)=a12(:,:)+(zp(1))**(k-1)*ddw3(:,k,:)
  fr1(:,:)=fr1(:,:)+(zp(1))**(k-1)*ddw1(:,k,:,1)
  fc1(:,:)=fc1(:,:)+(zp(1))**(k-1)*ddw1(:,k,:,2)
enddo
#endif


! upper boundary: z=+1
#if bc_up==0
! no-slip at z=+1
a21(:,:)=0.5d0*dw2(:,1,:)
a22(:,:)=0.5d0*dw3(:,1,:)
fr2(:,:)=0.5d0*dw1(:,1,:,1)
fc2(:,:)=0.5d0*dw1(:,1,:,2)
do k=2,nz
  a21(:,:)=a21(:,:)+dw2(:,k,:)
  a22(:,:)=a22(:,:)+dw3(:,k,:)
  fr2(:,:)=fr2(:,:)+dw1(:,k,:,1)
  fc2(:,:)=fc2(:,:)+dw1(:,k,:,2)
enddo
#else
! free-slip at z=+1
call dz(dw1,ddw1)
call dz_red(dw2,ddw2)
call dz_red(dw3,ddw3)
a21(:,:)=0.5d0*ddw2(:,1,:)
a22(:,:)=0.5d0*ddw3(:,1,:)
fr2(:,:)=0.5d0*ddw1(:,1,:,1)
fc2(:,:)=0.5d0*ddw1(:,1,:,2)
do k=2,nz
  a21(:,:)=a21(:,:)+ddw2(:,k,:)
  a22(:,:)=a22(:,:)+ddw3(:,k,:)
  fr2(:,:)=fr2(:,:)+ddw1(:,k,:,1)
  fc2(:,:)=fc2(:,:)+ddw1(:,k,:,2)
enddo
#endif


deallocate(dw1)
deallocate(dw2)
deallocate(dw3)
deallocate(ddw1)
deallocate(ddw2)
deallocate(ddw3)

! calculate complex valued coefficients A,B
allocate(A(spx,spy,2))
allocate(B(spx,spy,2))

do j=1,spy
  do i=1,spx
    det=(a11(i,j)*a22(i,j)-a21(i,j)*a12(i,j))
    A(i,j,1)=(-a22(i,j)*fr1(i,j)+a12(i,j)*fr2(i,j))/det
    A(i,j,2)=(-a22(i,j)*fc1(i,j)+a12(i,j)*fc2(i,j))/det
    B(i,j,1)=(a21(i,j)*fr1(i,j)-a11(i,j)*fr2(i,j))/det
    B(i,j,2)=(a21(i,j)*fc1(i,j)-a11(i,j)*fc2(i,j))/det
  enddo
enddo


deallocate(a11)
deallocate(a12)
deallocate(a21)
deallocate(a22)
deallocate(fr1)
deallocate(fr2)
deallocate(fc1)
deallocate(fc2)


! calculate w
do k=1,nz
  do j=1,spy
    do i=1,spx
      wc(i,k,j,1)=h(i,k,j,1)+A(i,j,1)*wa2(i,k,j)+B(i,j,1)*wa3(i,k,j)
      wc(i,k,j,2)=h(i,k,j,2)+A(i,j,2)*wa2(i,k,j)+B(i,j,2)*wa3(i,k,j)
    enddo
  enddo
enddo

deallocate(A)
deallocate(B)


! set zero mode for w (k^2=0)
if(rank.eq.0)then
  wc(1,:,1,1)=0.0d0
  wc(1,:,1,2)=0.0d0
endif


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_uv(omega,h1,h2)
! calculate u,v from the definition of vorticity and the
! continuity equation (w, omega_z already known at this step)

use commondata
use par_size
use wavenumber
use velocity
use sim_par

double precision, dimension(spx,nz,spy,2) :: omega,h1,h2,dw
double precision, allocatable, dimension(:) :: tempu,tempv
double precision :: beta

integer :: i,j,k

call dz(wc,dw)

do j=1,spy
  do k=1,nz
    do i=1,spx
      uc(i,k,j,1)=(-ky(j+cstart(3))*omega(i,k,j,2)-kx(i+cstart(1))*dw(i,k,j,2))/k2(i+cstart(1),j+cstart(3))
      uc(i,k,j,2)=(ky(j+cstart(3))*omega(i,k,j,1)+kx(i+cstart(1))*dw(i,k,j,1))/k2(i+cstart(1),j+cstart(3))
      vc(i,k,j,1)=(kx(i+cstart(1))*omega(i,k,j,2)-ky(j+cstart(3))*dw(i,k,j,2))/k2(i+cstart(1),j+cstart(3))
      vc(i,k,j,2)=(-kx(i+cstart(1))*omega(i,k,j,1)+ky(j+cstart(3))*dw(i,k,j,1))/k2(i+cstart(1),j+cstart(3))
    enddo
  enddo
enddo


! set value for k^2=0
if(rank.eq.0)then
  allocate(tempu(nz))
  allocate(tempv(nz))

  tempu(:)=-h1(1,:,1,1)/gamma
  tempv(:)=-h2(1,:,1,1)/gamma
  beta=1.0d0/gamma

  call helmholtz_rred(tempu,beta,p_u,q_u,r_u,zp)
  call helmholtz_rred(tempv,beta,p_v,q_v,r_v,zp)

  uc(1,:,1,1)=tempu(:)
  vc(1,:,1,1)=tempv(:)

  deallocate(tempu)
  deallocate(tempv)

  ! set imaginary part to 0
  uc(1,:,1,2) = 0.d0
  vc(1,:,1,2) = 0.d0
endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_omega(h1,h2,omega)

use commondata
use par_size
use sim_par
use velocity
use wavenumber

double precision, dimension(spx,nz,spy,2) :: h1,h2,omega
double precision, dimension(spx,spy) :: beta2
integer :: i,j


! set up RHS for Helmholtz equation and store it in omega
do i=1,spx
  omega(i,:,:,1)=-kx(i+cstart(1))*h2(i,:,:,2)
  omega(i,:,:,2)=kx(i+cstart(1))*h2(i,:,:,1)
enddo
do j=1,spy
  omega(:,:,j,1)=omega(:,:,j,1)+ky(j+cstart(3))*h1(:,:,j,2)
  omega(:,:,j,2)=omega(:,:,j,2)-ky(j+cstart(3))*h1(:,:,j,1)
enddo

omega=-omega/gamma


! set coefficient for LHS omega_z
do j=1,spy
  do i=1,spx
    beta2(i,j)=(1.0d0+gamma*k2(i+cstart(1),j+cstart(3)))/gamma
  enddo
enddo


! call Helmholtz solver
! omega_z will be saved in omega in output, omega is the RHS of the equation in input
call helmholtz(omega,beta2,p_o,q_o,r_o,zp)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_phi(hphi)

use commondata
use par_size
use phase_field
use sim_par
use wavenumber

double precision :: hphi(spx,nz,spy,2),beta2(spx,spy)

integer :: i,j

hphi=hphi*pe/(dt*ch**2)

do j=1,spy
  do i=1,spx
    beta2(i,j)=k2(i+cstart(1),j+cstart(3))+s_coeff/(2.0d0*ch**2)
  enddo
enddo

! solve for auxiliary variable
call helmholtz(hphi,beta2,[0.0d0,0.0d0],[1.0d0,1.0d0],[0.0d0,0.0d0],zp)



! solve for phi
call helmholtz(hphi,beta2,[0.0d0,0.0d0],[1.0d0,1.0d0],[0.0d0,0.0d0],zp)

phic=hphi

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_psi(hpsi)

use commondata
use par_size
use surfactant
use sim_par
use wavenumber

double precision :: hpsi(spx,nz,spy,2),beta2(spx,spy),gammapsi

integer :: i,j

gammapsi=dt*P_i/pe_psi

hpsi=-hpsi/gammapsi

do j=1,spy
  do i=1,spx
    beta2(i,j)=1.0d0/gammapsi+k2(i+cstart(1),j+cstart(3))
  enddo
enddo


call helmholtz(hpsi,beta2,[0.0d0,0.0d0],[1.0d0,1.0d0],[0.0d0,0.0d0],zp)

psic=hpsi


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_theta(htheta)

use commondata
use par_size
use sim_par
use wavenumber
use temperature

double precision :: htheta(spx,nz,spy,2),tmp(nz)
double precision :: gammatheta,beta2(spx,spy)

integer :: i,j

gammatheta=Re*Pr/dt

do j=1,spy
  do i=1,spx
    beta2(i,j)=k2(i+cstart(1),j+cstart(3))+gammatheta
  enddo
enddo

tmp=htheta(1,:,1,1)

call helmholtz(htheta,beta2,p_theta,q_theta,[0.0d0,0.0d0],zp)

thetac=htheta

if(rank.eq.0)then
 call helmholtz_rred(tmp,beta2(1,1),p_theta,q_theta,r_theta,zp)
 thetac(1,:,1,1)=tmp
endif



return
end
