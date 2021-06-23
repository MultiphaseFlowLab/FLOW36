subroutine get_center(s_drop,nstep,drop_count)

use commondata
use grid
use sim_par
use velocity

double precision :: sum_xy(nz)!,inertia_t(6)
double precision :: int_1,dx,dy,dz
double precision :: xg,yg,zg!,A ! que sont A Mx, My et Mz ?
double precision :: ug,vg,wg
!double precision, dimension(nx,nz,ny) :: phix,phiy,phiz

integer :: s_drop(nx,nz,ny)
integer :: nstep,drop_count
integer :: i,j,k,xboundary,yboundary,it,jt

character (LEN= 8):: T
write(T, '(i8.8)') nstep

dx=xl/dble(nx-1)
dy=yl/dble(ny-1)
! dx, dy are uniform in space, dz=dz(z), integral of phi=+1 on V can be reduced to
! a sum of integrals over z
sum_xy=sum(sum(dble(s_drop),3),1)

int_1=0.0d0; ! int_1 est le volume total de la droplet, qui est actualisé peu à peu
do k=1,nz-1
  int_1=int_1+(sum_xy(k)+sum_xy(k+1))*(z(k)-z(k+1))*0.5d0
enddo
int_1=int_1*dx*dy

!
xboundary=0
yboundary=0
! On regarde si la droplet et sur la frontière en x et / ou en y
! Si oui, on actualise xboundary et yboundary
do j=1,ny
  do k=1,nz
    do i=1,nx
      if(s_drop(i,k,j).eq.1)then
       if((i.eq.nx).or.(i.eq.1)) xboundary=1
       if((j.eq.ny).or.(j.eq.1)) yboundary=1
      endif
    enddo
  enddo
enddo

xg=0.0d0
yg=0.0d0
zg=0.0d0

ug=0.0d0
vg=0.0d0
wg=0.0d0

do j=1,ny
  do k=1,nz-1
    dz=z(k)-z(k+1)
    do i=1,nx
      if(xboundary.eq.1) then ! si sur la frontière en x, on fait une translation sur les indices par rapport à x
       it=mod(i-floor(real(nx)/4.0)+nx,nx)+1
      else
       it=i ! pas sur la frontière en x, on ne change pas les indices
      endif
      if(yboundary.eq.1) then ! si sur la frontière en y, on fait une translation sur les indices par rapport à y
       jt=mod(j-floor(real(ny)/4.0)+ny,ny)+1
      else
       jt=j ! pas sur la frontière en y, on ne change pas les indices
      endif

!s_drop(i,k,j) = 0 si pas dans la drop, = 1 si dans la drop
      xg=xg+dble(s_drop(i,k,j))*dx*dy*dz*x(it) ! it indice apres translation
      yg=yg+dble(s_drop(i,k,j))*dx*dy*dz*y(jt)
      zg=zg+dble(s_drop(i,k,j))*dx*dy*dz*z(k) !pas de translation selon z, indice reste le même

      ug=ug+dble(s_drop(i,k,j))*dble(u(i,k,j))*dx*dy*dz
      vg=vg+dble(s_drop(i,k,j))*dble(v(i,k,j))*dx*dy*dz
      wg=wg+dble(s_drop(i,k,j))*dble(w(i,k,j))*dx*dy*dz

    enddo
  enddo
enddo

xg=xg/int_1
yg=yg/int_1
zg=zg/int_1

ug=ug/int_1
vg=vg/int_1
wg=wg/int_1



! si on a translaté selon x
! on reranslate le centre de gravité
if(xboundary.eq.1)then
 if(floor(real(nx)/4.0).ge.1)then
  xg=xg+x(floor(real(nx)/4.0))
  if(xg.gt.xl) xg=xg-xl
 endif
endif


! si on a translaté selon y
! on reranslate le centre de gravité

if(yboundary.eq.1)then
  if(floor(real(ny)/4.0).ge.1)then
   yg=yg+y(floor(real(ny)/4.0))
   if(yg.gt.yl) yg=yg-yl
  endif
endif

!call normal(phix,phiy,phiz)

!call get_area(s_drop,A,phix,phiy,phiz,xg,yg,zg,Mx,My,Mz,xboundary,yboundary)

!call princ_axis(s_drop,xg,yg,zg,inertia_t,nstep,drop_count,xboundary,yboundary)

! comparison with volume of a sphere
! write(*,'(i8,5(f12.6))') nstep,xg,yg,zg,int_1,4.0d0/3.0d0*pi*0.3d0**3




! a decommenter!
open(70,file='./output/mass_center_global.dat',form='formatted',status='old',access='append')
! output in outer units
 ! write(66,'(i8,i6,5(es16.6))') nstep,drop_count,xg,yg,zg,int_1,A
! output in wall units
 ! à decommenter!
write(70,'(i8,i6,6x,4(es16.6))') nstep,drop_count,xg*Re,yg*Re,(zg+1.0d0)*Re, (int_1)*(Re**3)
!à décommenter !
close(70,status='keep')


! Différents fichiers (un pour chaque timestep)
! Chaque fichier contient pour le dit timestep la position du centre de gravité
! de chaque droplet présentes à ce timestep
!open(73,file='./output/mass_center_timestep_'//T//'.dat',form='formatted',status='old',access='append')
!write(73,'(i8,3(es16.6))') drop_count, xg*Re,yg*Re,(zg+1.0d0)*Re
!!!!write(73,'(i8,3(es16.6))') drop_count,xg, yg, zg
!close(73,status='keep')



!open(74,file='./output/velocity_center_timestep_'//T//'.dat',form='formatted',status='old',access='append')
!!!write(74,'(i8,3(es16.6))') drop_count, xg*Re,yg*Re,(zg+1.0d0)*Re
!write(74,'(i8,3(es16.6))') drop_count,ug, vg, wg
!close(74,status='keep')

open(74,file='./output/center_timestep_'//T//'.dat',form='formatted',status='old',access='append')
!write(74,'(i8,3(es16.6))') drop_count, xg*Re,yg*Re,(zg+1.0d0)*Re
write(74,'(i8,7(es16.6))') drop_count, xg*Re, yg*Re, (zg+1.0d0)*Re, ug, vg, wg, (int_1)*(Re**3)
close(74,status='keep')


! a decommenter!
!open(66,file='./output/mass_center.dat',form='formatted',status='old',access='append')
! output in outer units
 ! write(66,'(i8,i6,5(es16.6))') nstep,drop_count,xg,yg,zg,int_1,A
! output in wall units
 ! à decommenter!
!write(66,'(i8,i6,6x,3(es16.6))') nstep,drop_count,xg*Re,yg*Re,(zg+1.0d0)*Re
!à décommenter !
!close(66,status='keep')




!open(66,file='./output/mass_center.dat',form='formatted',status='old',access='append')
!! output in outer units
 !! write(66,'(i8,i6,5(es16.6))') nstep,drop_count,xg,yg,zg,int_1,A
!! output in wall units
! write(66,'(i8,i6,8(es16.6))') nstep,drop_count,xg*Re,yg*Re,(zg+1.0d0)*Re,int_1*Re**3,A*Re**2,Mx*Re**3,My*Re**3,Mz*Re**3
!close(66,status='keep')

!open(67,file='./output/inertia_tensor.dat',form='formatted',status='old',access='append')
 !write(67,'(2(i6),6(es16.6))') nstep,drop_count,inertia_t(1:6)
!close(67,status='keep')

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
