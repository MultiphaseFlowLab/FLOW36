subroutine dist_function(nstep)

use commondata

integer :: nstep
integer :: i,j,k,m,ip,jp,kp,pos(1)
double precision :: step,phi0,phip,xp,yp,zp
double precision, dimension(nx,nz,ny) :: dist
double precision, dimension(nx,nz,ny) :: modg,derx,dery,derz


! explore neighbouring volume
! call explore(dist)

! ! follow normal at each point until interface is found
! call mod_grad(phi,modg,derx,dery,derz)
!
! step=(y(2)-y(1))/10.0d0
! dist=0.0d0
!
! do j=1,ny
!  do k=1,nz
!   do i=1,nx
!    phi0=phi(i,k,j)
!    phip=phi0
!    m=0
!    do while((dble(m)*step.lt.5.0d0).or.(phip*phi0.le.0.0d0)) ! loop until interface is found or interface is farther than ...
!     m=m+1 ! m=0 already considered, if phi0=0 loop is skipped
!     ! direction of investigation depends on local value of phi
!     xp=x(i)-phi0/abs(phi0)*dble(m)*step*derx(i,k,j)
!     yp=y(j)-phi0/abs(phi0)*dble(m)*step*dery(i,k,j)
!     zp=z(k)-phi0/abs(phi0)*dble(m)*step*derz(i,k,j)
!     ! consider periodicity and return if outside z
!     if(xp.gt.lx)then
!      xp=xp-lx
!     elseif(xp.lt.0)then
!      xp=xp+lx
!     endif
!     if(yp.gt.ly)then
!      yp=yp-ly
!     elseif(yp.lt.0)then
!      yp=yp+ly
!     endif
!     if((zp.gt.lz).or.(zp.lt.0))then
!      ! put some condition on distance
!      exit
!     endif
!
!     pos=minloc(dabs(xp-x))
!     ip=pos(1)
!     if(xp.lt.x(ip)) ip=ip-1
!     pos=minloc(dabs(yp-y))
!     jp=pos(1)
!     if(yp.lt.y(jp)) jp=jp-1
!     pos=minloc(dabs(zp-z))
!     kp=pos(1)
!     if(kp.lt.z(kp)) kp=kp+1
!
!     call trilinear(phi,ip,jp,kp,xp,yp,zp,phip)
!    enddo
!    dist(i,k,j)=dble(m)*step
!   enddo
!  enddo
! enddo


! calculate distance knowing that profile is tanh
do j=1,ny
 do k=1,nz
  do i=1,nx
   dist(i,k,j)=dsqrt(2.0d0)*ch*atanh(phi(i,k,j))
  enddo
 enddo
enddo
! output in minus units
! problem: oscillations in phi profile, either solve d phi/dt=f_profile or cut fluctuations
! solve d phi/dt=f_profile better solution, might be costly


call generate_output(nstep,dist)

! redistancing equation unstable AF
! call mod_grad(phi,modg,derx,dery,derz)
! control=maxval(abs(modg))
!
! do while(abs(control-1.0d0).gt.0.1)
!   ! already update phi value
!   phi=phi-tt*sgn*(modg-1.0d0)
!   ! call sign_f(phi,sgn)  !!!! keep original interface position
!   call mod_grad(phi,modg,derx,dery,derz)
!   control=maxval(abs(modg))
!   write(*,*) control
!   ! if(control.lt.5.0d0) tt=tt*2
! enddo



return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mod_grad(fun,modg,derx,dery,derz)

use commondata
use wavenumber

double precision, dimension(nx,nz,ny) :: modg,derx,dery,derz,fun
double precision, dimension(nx/2+1,nz,ny,2) :: tmpc,func

integer :: i

call phys_to_spectral(fun,func,1)

do i=1,nx/2+1
  tmpc(i,:,:,1)=-kx(i)*func(i,:,:,2)
  tmpc(i,:,:,2)=+kx(i)*func(i,:,:,1)
enddo
call spectral_to_phys(tmpc,derx,1)

do i=1,ny
  tmpc(:,:,i,1)=-ky(i)*func(:,:,i,2)
  tmpc(:,:,i,2)=+ky(i)*func(:,:,i,1)
enddo
call spectral_to_phys(tmpc,dery,1)

call dz(func,tmpc)
call spectral_to_phys(tmpc,derz,1)

modg=(derx**2+dery**2+derz**2)**0.5d0

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sign_f(fun,sgn)

use commondata

double precision, dimension(nx,nz,ny) :: sgn,fun

integer :: i,k,j

do j=1,ny
 do k=1,nz
  do i=1,nx
    if(fun(i,k,j).gt.0.0d0)then
      sgn(i,k,j)=+1.0d0
    elseif(fun(i,k,j).eq.0.0d0)then
      sgn(i,k,j)=0.0d0
    elseif(fun(i,k,j).lt.0.0d0)then
      sgn(i,k,j)=-1.0d0
    endif
  enddo
 enddo
enddo

! sgn=fun/dabs(fun)

! write(*,*) maxval(sgn),minval(sgn)

return
end
