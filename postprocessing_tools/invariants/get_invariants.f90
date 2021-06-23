subroutine get_invariants(nstep)

use commondata
use wavenumber

integer :: nstep
integer :: i,k,j

double precision, dimension(nx,nz,ny) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dupdz,dvpdz,dwpdz
double precision, dimension(nz) :: um,vm,wm
double precision, allocatable, dimension(:,:,:,:) :: tmpcu,tmpcv,tmpcw
double precision :: Q_v,R_v,Q_s,R_s,Q_omega

character(len=8) :: step

write(step,'(i8.8)') nstep

allocate(tmpcu(nx/2+1,nz,ny,2))
allocate(tmpcv(nx/2+1,nz,ny,2))
allocate(tmpcw(nx/2+1,nz,ny,2))


! uc,vc,wc already available in spectral space, no need for transforms
! calculation of velocity gradient tensor
do i=1,nx/2+1
  tmpcu(i,:,:,1)=-kx(i)*uc(i,:,:,2)
  tmpcu(i,:,:,2)=+kx(i)*uc(i,:,:,1)
  tmpcv(i,:,:,1)=-kx(i)*vc(i,:,:,2)
  tmpcv(i,:,:,2)=+kx(i)*vc(i,:,:,1)
  tmpcw(i,:,:,1)=-kx(i)*wc(i,:,:,2)
  tmpcw(i,:,:,2)=+kx(i)*wc(i,:,:,1)
enddo
call spectral_to_phys(tmpcu,dudx,0)
call spectral_to_phys(tmpcv,dvdx,0)
call spectral_to_phys(tmpcw,dwdx,0)

do j=1,ny
  tmpcu(:,:,j,1)=-ky(j)*uc(:,:,j,2)
  tmpcu(:,:,j,2)=+ky(j)*uc(:,:,j,1)
  tmpcv(:,:,j,1)=-ky(j)*vc(:,:,j,2)
  tmpcv(:,:,j,2)=+ky(j)*vc(:,:,j,1)
  tmpcw(:,:,j,1)=-ky(j)*wc(:,:,j,2)
  tmpcw(:,:,j,2)=+ky(j)*wc(:,:,j,1)
enddo
call spectral_to_phys(tmpcu,dudy,0)
call spectral_to_phys(tmpcv,dvdy,0)
call spectral_to_phys(tmpcw,dwdy,0)

call dz(uc,tmpcu)
call dz(vc,tmpcw)
call dz(wc,tmpcw)
call spectral_to_phys(tmpcu,dudz,0)
call spectral_to_phys(tmpcv,dvdz,0)
call spectral_to_phys(tmpcw,dwdz,0)


! derivatives of mean velocities in homogenous directions are zero, only wall-normal derivatives are calculated
um=sum(sum(u,3),1)/dble(nx*ny)
vm=sum(sum(v,3),1)/dble(nx*ny)
wm=sum(sum(w,3),1)/dble(nx*ny)

! u,v,w,uc,vc,wc now contain fluctuating component!!!
do k=1,nz
  u(:,k,:)=u(:,k,:)-um(k)
  v(:,k,:)=v(:,k,:)-vm(k)
  w(:,k,:)=w(:,k,:)-wm(k)
enddo

call phys_to_spectral(u,uc,0)
call phys_to_spectral(v,vc,0)
call phys_to_spectral(w,wc,0)
call dz(uc,tmpcu)
call dz(vc,tmpcw)
call dz(wc,tmpcw)
call spectral_to_phys(tmpcu,dupdz,0)
call spectral_to_phys(tmpcv,dvpdz,0)
call spectral_to_phys(tmpcw,dwpdz,0)
! velocity gradient tensor calculated + fluctuating components

deallocate(tmpcu,tmpcv,tmpcw)

! open(666,file='./output/inv_in_'//step//'.dat',form='formatted',status='new')
! open(667,file='./output/inv_out_'//step//'.dat',form='formatted',status='new')
! write(666,'(6(a18))') 'Q_v','R_v','Q_s','R_s','Q_omega','inside drop'
! write(667,'(6(a18))') 'Q_v','R_v','Q_s','R_s','Q_omega','outside drop'

do j=1,ny
 do k=1,nz
  do i=1,nx
   ! Q_v
   Q_v=-(-dudx(i,k,j)*dvdy(i,k,j)-dudx(i,k,j)*dwdz(i,k,j)-dvdy(i,k,j)*dwdz(i,k,j) &
 &                     +dvdx(i,k,j)*dudy(i,k,j)+dudz(i,k,j)*dwdx(i,k,j)+dvdz(i,k,j)*dwdy(i,k,j))
   ! R_v
   R_v=-((dudx(i,k,j))**3+(dvdy(i,k,j))**3+(dwdz(i,k,j))**3+3.0d0*( &
 &                     dudx(i,k,j)*dudy(i,k,j)*dvdx(i,k,j)+dudx(i,k,j)*dudz(i,k,j)*dwdx(i,k,j)+ &
 &                     dudy(i,k,j)*dvdy(i,k,j)*dvdx(i,k,j)+dudy(i,k,j)*dvdz(i,k,j)*dwdx(i,k,j)+ &
 &                     dudz(i,k,j)*dwdy(i,k,j)*dvdx(i,k,j)+dudz(i,k,j)*dwdz(i,k,j)*dwdx(i,k,j)+ &
 &                     dvdy(i,k,j)*dvdz(i,k,j)*dwdy(i,k,j)+dvdz(i,k,j)*dwdz(i,k,j)*dwdy(i,k,j)))/3.0d0
   ! Q_s
   Q_s=-((2.0d0*dudx(i,k,j))**2+2.0d0*(dudy(i,k,j)+dvdx(i,k,j))**2+2.0d0*(dudz(i,k,j)+dwdx(i,k,j))**2 &
 &                      +(2.0d0*dvdy(i,k,j))**2+2.0d0*(dvdz(i,k,j)+dwdy(i,k,j))**2 &
 &                      +(2.0d0*dwdz(i,k,j))**2)/8.0d0
   ! R_s
   R_s=-((2.0d0*dudx(i,k,j))**3+(2.0d0*dvdy(i,k,j))**3+(2.0d0*dwdz(i,k,j))**3+3.0d0*( &
 &          dudx(i,k,j)*(dudy(i,k,j)+dvdx(i,k,j))**2+dudx(i,k,j)*(dudz(i,k,j)+dwdx(i,k,j))**2+ &
 &          dvdy(i,k,j)*(dudy(i,k,j)+dvdx(i,k,j))**2+dwdz(i,k,j)*(dudz(i,k,j)+dwdx(i,k,j))**2+ &
 &          dvdy(i,k,j)*(dvdz(i,k,j)+dwdy(i,k,j))**2+dwdz(i,k,j)*(dvdz(i,k,j)+dwdy(i,k,j))**2+ &
 &          2.0d0*(dvdx(i,k,j)+dudy(i,k,j))*(dwdy(i,k,j)+dvdz(i,k,j))*(dwdx(i,k,j)+dudz(i,k,j))))/(8.0d0*3.0d0)
   ! Q_omega
   Q_omega=2.0d0*((dudy(i,k,j)-dvdx(i,k,j))**2+(dudz(i,k,j)-dwdx(i,k,j))**2+(dwdy(i,k,j)-dvdz(i,k,j))**2)/8.0d0

   ! if(phi(i,k,j).ge.0.0d0)then
   !  write(666,'(5(es18.6))') Q_v,R_v,Q_s,R_s,Q_omega
   !  samples_in=samples_in+1
   ! else
   !  write(667,'(5(es18.6))') Q_v,R_v,Q_s,R_s,Q_omega
   !  samples_out=samples_out+1
   ! endif
   call compute_jpdf(Q_v,R_v,Q_s,R_s,Q_omega,phi(i,k,j))
  enddo
 enddo
enddo

! close(666,status='keep')
! close(667,status='keep')


return
end
