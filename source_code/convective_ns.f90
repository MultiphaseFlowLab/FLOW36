subroutine convective_ns(s1,s2,s3)

use commondata
use velocity
use par_size
use wavenumber
use phase_field

double precision :: s1(spx,nz,spy,2), s2(spx,nz,spy,2), s3(spx,nz,spy,2)
double precision, allocatable, dimension(:,:,:) :: uu,uv,uw,vv,vw,ww
double precision, allocatable, dimension(:,:,:,:) :: uuc,uvc,uwc,vvc,vwc,wwc
double precision, allocatable, dimension(:,:,:,:) :: uucx,uvcy,uwcz, uvcx,vvcy,vwcz, uwcx,vwcy,wwcz
double precision :: phif


integer :: indx,indy
integer :: i,j,k

! transform variables back to physical space and perform dealiasing
call spectral_to_phys(uc,u,1,1)
call spectral_to_phys(vc,v,1,1)
call spectral_to_phys(wc,w,1,1)

call courant_check(u,v,w)

allocate(uu(nx,fpz,fpy))
allocate(uv(nx,fpz,fpy))
allocate(uw(nx,fpz,fpy))
allocate(vv(nx,fpz,fpy))
allocate(vw(nx,fpz,fpy))
allocate(ww(nx,fpz,fpy))

allocate(uuc(spx,nz,spy,2))
allocate(uvc(spx,nz,spy,2))
allocate(uwc(spx,nz,spy,2))
allocate(vvc(spx,nz,spy,2))
allocate(vwc(spx,nz,spy,2))
allocate(wwc(spx,nz,spy,2))


! form products uu,uv,uw, ...
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      uu(i,k,j)=u(i,k,j)*u(i,k,j)
      uv(i,k,j)=u(i,k,j)*v(i,k,j)
      uw(i,k,j)=u(i,k,j)*w(i,k,j)
      vv(i,k,j)=v(i,k,j)*v(i,k,j)
      vw(i,k,j)=v(i,k,j)*w(i,k,j)
      ww(i,k,j)=w(i,k,j)*w(i,k,j)
    enddo
  enddo
enddo


! transform products uu,uv, ...  in spectral space
call phys_to_spectral(uu,uuc,1,1)
call phys_to_spectral(uv,uvc,1,1)
call phys_to_spectral(uw,uwc,1,1)
call phys_to_spectral(vv,vvc,1,1)
call phys_to_spectral(vw,vwc,1,1)
call phys_to_spectral(ww,wwc,1,1)

deallocate(uu)
deallocate(uv)
deallocate(uw)
deallocate(vv)
deallocate(vw)
deallocate(ww)


! take derivatives
allocate(uucx(spx,nz,spy,2))
allocate(uvcy(spx,nz,spy,2))
allocate(uwcz(spx,nz,spy,2))
allocate(uvcx(spx,nz,spy,2))
allocate(vvcy(spx,nz,spy,2))
allocate(vwcz(spx,nz,spy,2))
allocate(uwcx(spx,nz,spy,2))
allocate(vwcy(spx,nz,spy,2))
allocate(wwcz(spx,nz,spy,2))

! x derivatives
indx=cstart(1)
do i=1,spx
  uucx(i,:,:,1)=-uuc(i,:,:,2)*kx(indx+i)
  uucx(i,:,:,2)=uuc(i,:,:,1)*kx(indx+i)
  uvcx(i,:,:,1)=-uvc(i,:,:,2)*kx(indx+i)
  uvcx(i,:,:,2)=uvc(i,:,:,1)*kx(indx+i)
  uwcx(i,:,:,1)=-uwc(i,:,:,2)*kx(indx+i)
  uwcx(i,:,:,2)=uwc(i,:,:,1)*kx(indx+i)
enddo


! y derivatives
indy=cstart(3)
do j=1,spy
  uvcy(:,:,j,1)=-uvc(:,:,j,2)*ky(indy+j)
  uvcy(:,:,j,2)=uvc(:,:,j,1)*ky(indy+j)
  vvcy(:,:,j,1)=-vvc(:,:,j,2)*ky(indy+j)
  vvcy(:,:,j,2)=vvc(:,:,j,1)*ky(indy+j)
  vwcy(:,:,j,1)=-vwc(:,:,j,2)*ky(indy+j)
  vwcy(:,:,j,2)=vwc(:,:,j,1)*ky(indy+j)
enddo


! z derivatives
call dz(uwc,uwcz)
call dz(vwc,vwcz)
call dz(wwc,wwcz)



deallocate(uuc)
deallocate(uvc)
deallocate(uwc)
deallocate(vvc)
deallocate(vwc)
deallocate(wwc)


! form non linear terms
s1=-(uucx+uvcy+uwcz)
s2=-(uvcx+vvcy+vwcz)
s3=-(uwcx+vwcy+wwcz)


deallocate(uucx)
deallocate(uvcy)
deallocate(uwcz)
deallocate(uvcx)
deallocate(vvcy)
deallocate(vwcz)
deallocate(uwcx)
deallocate(vwcy)
deallocate(wwcz)


! conditional compilation, add phase variable contribution for non-matched densities
#define match_dens matched_density
#define phiflag phicompflag
#if phiflag == 1
#if match_dens != 1
! case for non-matched densities

allocate(uu(nx,fpz,fpy))
allocate(vv(nx,fpz,fpy))
allocate(ww(nx,fpz,fpy))

call spectral_to_phys(s1,uu,1,0)
call spectral_to_phys(s2,vv,1,0)
call spectral_to_phys(s3,ww,1,0)

call spectral_to_phys(phic,phi,1,0)

#if match_dens == 2
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      uu(i,k,j)=(1.0d0+0.5d0*(1.0d0/rhor-1.0d0)*(1.0d0-phif))*uu(i,k,j)
      vv(i,k,j)=(1.0d0+0.5d0*(1.0d0/rhor-1.0d0)*(1.0d0-phif))*vv(i,k,j)
      ww(i,k,j)=(1.0d0+0.5d0*(1.0d0/rhor-1.0d0)*(1.0d0-phif))*ww(i,k,j)
    enddo
  enddo
enddo
#else
do j=1,fpy
  do k=1,fpz
    do i=1,nx
! phi filtered to remove overshooting when calculating rho
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      uu(i,k,j)=(1.0d0+0.5d0*(rhor-1.0d0)*(phif+1.0d0))*uu(i,k,j)
      vv(i,k,j)=(1.0d0+0.5d0*(rhor-1.0d0)*(phif+1.0d0))*vv(i,k,j)
      ww(i,k,j)=(1.0d0+0.5d0*(rhor-1.0d0)*(phif+1.0d0))*ww(i,k,j)
    enddo
  enddo
enddo
#endif

call phys_to_spectral(uu,s1,1,0)
call phys_to_spectral(vv,s2,1,0)
call phys_to_spectral(ww,s3,1,0)

deallocate(uu)
deallocate(vv)
deallocate(ww)

#endif
#endif

return
end
