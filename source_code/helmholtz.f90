subroutine helmholtz(f,beta2,p,q,r,z_p)
! solve equation
! y"-beta2*y=f
! with boundary conditions
! p(i)*y(z(i))+q(i)*y'(z(i))=r(i)       @ z(i) , i=1,2
!
! Same algorithm and operations as the previous version (system assembled in the 3D arrays
! a,b,c,d,e,h and solved by gauss_solver), fused in one kernel per (i,j) column:
! - a(k-2)=-beta2*k, c(k-2)=-beta2*(k-2) and the initial b(k-2) are computed on the fly
! - d,e depend only on k: 1D arrays d0,e0
! - during the elimination (k=nz..3) the modified b(k-4), d(k-2), e(k-2), h(k-2) are used two
!   steps later (index of the same parity): they are kept in parity slots bq,dq,eq,hq
! - only the final h(k) and b(k-2) needed by the forward substitution are stored (hs,bs)
! - the solution is written directly in f (no f=h copy)

use commondata
use par_size

double precision, dimension(spx,nz,spy,2) :: f
double precision, dimension(spx,spy) :: beta2
double precision, dimension(2) :: p,q,r,z_p
! work arrays, allocated once (managed memory: no allocation or transfer at every call)
double precision, allocatable, save :: hs(:,:,:,:),bs(:,:,:),d0(:),e0(:)
double precision :: t,dert,bt,ak2,bk2,rc,rd,re,d1,d2,e1,e2
double precision :: bq(0:1),dq(0:1),eq(0:1),hq(0:1,2),xs(0:1,2),hk(2),h1(2),h2(2),x1(2),x2(2)
integer :: i,j,k,kk,pk,c

if(.not.allocated(hs))then
 allocate(hs(spx,nz,spy,2),bs(spx,nz-2,spy),d0(nz),e0(nz))
endif

! first and second row (boundary conditions), see Canuto et al. 2006, pag. 85
!$acc parallel loop private(t,dert)
do k=1,nz
 if(k.eq.1)then
  d0(1)=0.5d0*(p(1)*1.0d0+q(1)*0.0d0)
  e0(1)=0.5d0*(p(2)*1.0d0+q(2)*0.0d0)
 else
  t=dble((z_p(1))**(k-1))
  dert=dble((z_p(1))**(k)*(k-1)**2)
  d0(k)=p(1)*t+q(1)*dert
  t=dble((z_p(2))**(k-1))
  dert=dble((z_p(2))**(k)*(k-1)**2)
  e0(k)=p(2)*t+q(2)*dert
 endif
enddo

!$acc parallel loop collapse(2) private(bq,dq,eq,hq,xs,hk,h1,h2,x1,x2)
do j=1,spy
 do i=1,spx
  bt=beta2(i,j)
  ! slots for the first two steps (k=nz, nz-1): values not modified yet
  do kk=nz-1,nz
   pk=mod(kk,2)
   bq(pk)=4.0d0*dble(kk*(kk-1)*(kk-2))+dble(2*(kk-1))*bt
   dq(pk)=d0(kk)
   eq(pk)=e0(kk)
  enddo
  do c=1,2
   hq(mod(nz,2),c)=dble(nz)*f(i,nz-2,j,c)-dble(2*(nz-1))*f(i,nz,j,c)
   hq(mod(nz-1,2),c)=dble(nz-1)*f(i,nz-3,j,c)-dble(2*(nz-2))*f(i,nz-1,j,c)
   h1(c)=r(1)
   h2(c)=r(2)
  enddo
  ! elimination, k=nz..3 (gauss_solver: upper diagonal c only for k>=5)
  do k=nz,3,-1
   pk=mod(k,2)
   bk2=bq(pk)
   bs(i,k-2,j)=bk2
   ak2=-bt*dble(k)
   do c=1,2
    hk(c)=hq(pk,c)
    hs(i,k,j,c)=hk(c)
   enddo
   if(k.ge.5)then
    ! cancel upper diagonal c
    rc=(-bt*dble(k-4))/bk2
    kk=k-2
    bq(pk)=(4.0d0*dble(kk*(kk-1)*(kk-2))+dble(2*(kk-1))*bt)-rc*ak2
    do c=1,2
     hq(pk,c)=(dble(kk)*f(i,kk-2,j,c)-dble(2*(kk-1))*f(i,kk,j,c)+dble(kk-2)*f(i,kk+2,j,c))-rc*hk(c)
    enddo
   endif
   ! cancel part of row d
   rd=dq(pk)/bk2
   dq(pk)=d0(k-2)-rd*ak2
   ! cancel part of row e
   re=eq(pk)/bk2
   eq(pk)=e0(k-2)-re*ak2
   do c=1,2
    h1(c)=h1(c)-rd*hk(c)
    h2(c)=h2(c)-re*hk(c)
   enddo
  enddo
  ! remains only 2x2 array (rows d,e reduced to their first two entries)
  d1=dq(1)
  d2=dq(0)
  e1=eq(1)
  e2=eq(0)
  do c=1,2
   x1(c)=(h1(c)-h2(c)*d2/e2)/(d1-e1*d2/e2)
   x2(c)=(h2(c)-x1(c)*e1)/(e2)
   f(i,1,j,c)=x1(c)
   f(i,2,j,c)=x2(c)
   xs(1,c)=x1(c)
   xs(0,c)=x2(c)
  enddo
  ! forward-substitute variable from k=3,nz, solution stored in f
  do k=3,nz
   pk=mod(k,2)
   do c=1,2
    xs(pk,c)=(hs(i,k,j,c)-(-bt*dble(k))*xs(pk,c))/bs(i,k-2,j)
    f(i,k,j,c)=xs(pk,c)
   enddo
  enddo
 enddo
enddo

return
end

! shape of the array resulting from the solver, only saves the 2 rows and the 3 diagonals
! a1, ... , an : Chebyshev coefficients (unknowns) of y
! r1,r2,b3, ... , bn : Chebyshev coefficients of f
! _                                   _   _  _   _  _
!| d d d d d d d d d d d d d d d d d d | | a1 | | r1 |
!| e e e e e e e e e e e e e e e e e e | | a2 | | r2 |
!| a 0 b 0 c 0 0 0 0 0 0 0 0 0 0 0 0 0 | | a3 | | b3 |
!| 0 a 0 b 0 c 0 0 0 0 0 0 0 0 0 0 0 0 | | a4 | | b4 |
!| 0 0 a 0 b 0 c 0 0 0 0 0 0 0 0 0 0 0 | | a5 | | b5 |
!| ................................... |*| .. |=| .. |
!| ................................... | | .. | | .. |
!| 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b 0 c 0 | | .. | | .. |
!| 0 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b 0 c | | .. | | .. |
!| 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b 0 | | .. | | .. |
!|_0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b_| |_an_| |_bn_|
!

subroutine gauss_solver(a,b,c,d,e,f)

use commondata
use par_size

double precision, dimension(spx,nz,spy,2) :: f
double precision :: a(spx,nz-2,spy),b(spx,nz-2,spy),c(spx,nz-4,spy),d(spx,nz,spy),e(spx,nz,spy)
double precision :: rc,rd,re

integer :: i,j,k

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    do k=nz,5,-1
      ! cancel upper diagonal c
      rc=c(i,k-4,j)/b(i,k-2,j)
      b(i,k-4,j)=b(i,k-4,j)-rc*a(i,k-2,j)
      f(i,k-2,j,1)=f(i,k-2,j,1)-rc*f(i,k,j,1)
      f(i,k-2,j,2)=f(i,k-2,j,2)-rc*f(i,k,j,2)
      ! cancel part of row d
      rd=d(i,k,j)/b(i,k-2,j)
      d(i,k-2,j)=d(i,k-2,j)-rd*a(i,k-2,j)
      f(i,1,j,1)=f(i,1,j,1)-rd*f(i,k,j,1)
      f(i,1,j,2)=f(i,1,j,2)-rd*f(i,k,j,2)
      ! cancel part of row e
      re=e(i,k,j)/b(i,k-2,j)
      e(i,k-2,j)=e(i,k-2,j)-re*a(i,k-2,j)
      f(i,2,j,1)=f(i,2,j,1)-re*f(i,k,j,1)
      f(i,2,j,2)=f(i,2,j,2)-re*f(i,k,j,2)
    enddo
! same as before for arrays d,e,f, but here no more upper diagonal
    do k=4,3,-1
      ! cancel part of row d
      rd=d(i,k,j)/b(i,k-2,j)
      d(i,k-2,j)=d(i,k-2,j)-rd*a(i,k-2,j)
      f(i,1,j,1)=f(i,1,j,1)-rd*f(i,k,j,1)
      f(i,1,j,2)=f(i,1,j,2)-rd*f(i,k,j,2)
      ! cancel part of row e
      re=e(i,k,j)/b(i,k-2,j)
      e(i,k-2,j)=e(i,k-2,j)-re*a(i,k-2,j)
      f(i,2,j,1)=f(i,2,j,1)-re*f(i,k,j,1)
      f(i,2,j,2)=f(i,2,j,2)-re*f(i,k,j,2)
    enddo
    ! remains only 2x2 array, diagonal and lower diagonal
    !  |d(1), d(2)|  |f(1)|
    !  |e(1), e(2)|  |f(2)|
    ! variable stored in f(1)
    f(i,1,j,1)=(f(i,1,j,1)-f(i,2,j,1)*d(i,2,j)/e(i,2,j))/(d(i,1,j)-e(i,1,j)*d(i,2,j)/e(i,2,j))
    f(i,1,j,2)=(f(i,1,j,2)-f(i,2,j,2)*d(i,2,j)/e(i,2,j))/(d(i,1,j)-e(i,1,j)*d(i,2,j)/e(i,2,j))
    ! variable stored in f(2)
    f(i,2,j,1)=(f(i,2,j,1)-f(i,1,j,1)*e(i,1,j))/(e(i,2,j))
    f(i,2,j,2)=(f(i,2,j,2)-f(i,1,j,2)*e(i,1,j))/(e(i,2,j))

    ! forward-substitute variable from k=2,Nz
    ! solution stored in f
    do k=3,nz
      f(i,k,j,1)=(f(i,k,j,1)-a(i,k-2,j)*f(i,k-2,j,1))/b(i,k-2,j)
      f(i,k,j,2)=(f(i,k,j,2)-a(i,k-2,j)*f(i,k-2,j,2))/b(i,k-2,j)
    enddo

  enddo
enddo


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine helmholtz_fg(f,beta2,p,q,r,z_p)
! solve equation
! y"-beta2*y=f
! with boundary conditions
! p(i)*y(z(i))+q(i)*y'(z(i))=r(i)       @ z(i) , i=1,2

use commondata
use par_size
use dual_grid

double precision, dimension(spxpsi,npsiz,spypsi,2) :: f,h
double precision, dimension(spxpsi,spypsi) :: beta2
double precision, dimension(2) :: p,q,r,z_p
double precision, dimension(spxpsi,npsiz-2,spypsi) :: a,b
double precision, dimension(spxpsi,npsiz-4,spypsi) :: c
double precision, dimension(spxpsi,npsiz,spypsi) :: d,e
double precision :: t0,dt0,t,dert

integer :: i,j,k

! system has line 1 and 2 full and then it is tridiagonal from line 3 to Nz
! a is the lower diagonal, b the diagonal, c the upper diagonal, d the first row, e the second row
! in the diagonal are excluded the first 2 rows (already in the arrays d and e)


! assemble diagonals
!$acc parallel loop collapse(2)
do j=1,spypsi
  do k=3,npsiz
    do i=1,spxpsi
      a(i,k-2,j)=-beta2(i,j)*dble(k)
      b(i,k-2,j)=4.0d0*dble(k*(k-1)*(k-2))+dble(2*(k-1))*beta2(i,j)
    enddo
  enddo
enddo

do j=1,spypsi
  do k=3,npsiz-2
    do i=1,spxpsi
      c(i,k-2,j)=-beta2(i,j)*dble(k-2)
    enddo
  enddo
enddo


! assemble 1st row
t0=1.0d0
dt0=0.0d0
d(:,1,:)=0.5d0*(p(1)*t0+q(1)*dt0)
do k=2,npsiz
! see Canuto et al. 2006, pag. 85
  t=dble((z_p(1))**(k-1))
  dert=dble((z_p(1))**(k)*(k-1)**2)
  d(:,k,:)=p(1)*t+q(1)*dert
enddo


! assemble 2nd row
t0=1.0d0
dt0=0.0d0
e(:,1,:)=0.5d0*(p(2)*t0+q(2)*dt0)
do k=2,npsiz
! see Canuto et al. 2006, pag. 85
  t=dble((z_p(2))**(k-1))
  dert=dble((z_p(2))**(k)*(k-1)**2)
  e(:,k,:)=p(2)*t+q(2)*dert
enddo


! assemble RHS of equation
h(:,1,:,:)=r(1)
h(:,2,:,:)=r(2)

do k=3,npsiz-2
  h(:,k,:,:)=dble(k)*f(:,k-2,:,:)-dble(2*(k-1))*f(:,k,:,:)+dble(k-2)*f(:,k+2,:,:)
enddo
h(:,npsiz-1,:,:)=dble(npsiz-1)*f(:,npsiz-3,:,:)-dble(2*(npsiz-2))*f(:,npsiz-1,:,:)
h(:,npsiz,:,:)=dble(npsiz)*f(:,npsiz-2,:,:)-dble(2*(npsiz-1))*f(:,npsiz,:,:)


call gauss_solver_fg(a,b,c,d,e,h)

f=h

return
end

! shape of the array resulting from the solver, only saves the 2 rows and the 3 diagonals
! a1, ... , an : Chebyshev coefficients (unknowns) of y
! r1,r2,b3, ... , bn : Chebyshev coefficients of f
! _                                   _   _  _   _  _
!| d d d d d d d d d d d d d d d d d d | | a1 | | r1 |
!| e e e e e e e e e e e e e e e e e e | | a2 | | r2 |
!| a 0 b 0 c 0 0 0 0 0 0 0 0 0 0 0 0 0 | | a3 | | b3 |
!| 0 a 0 b 0 c 0 0 0 0 0 0 0 0 0 0 0 0 | | a4 | | b4 |
!| 0 0 a 0 b 0 c 0 0 0 0 0 0 0 0 0 0 0 | | a5 | | b5 |
!| ................................... |*| .. |=| .. |
!| ................................... | | .. | | .. |
!| 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b 0 c 0 | | .. | | .. |
!| 0 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b 0 c | | .. | | .. |
!| 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b 0 | | .. | | .. |
!|_0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a 0 b_| |_an_| |_bn_|
!

subroutine gauss_solver_fg(a,b,c,d,e,f)

use commondata
use par_size
use dual_grid

double precision, dimension(spxpsi,npsiz,spypsi,2) :: f
double precision, dimension(spxpsi,npsiz-2,spypsi) :: a,b
double precision, dimension(spxpsi,npsiz-4,spypsi) :: c
double precision, dimension(spxpsi,npsiz,spypsi) :: d,e
double precision :: rc,rd,re

integer :: i,j,k


do j=1,spypsi
  do i=1,spxpsi
    do k=npsiz,5,-1
      ! cancel upper diagonal c
      rc=c(i,k-4,j)/b(i,k-2,j)
      b(i,k-4,j)=b(i,k-4,j)-rc*a(i,k-2,j)
      f(i,k-2,j,1)=f(i,k-2,j,1)-rc*f(i,k,j,1)
      f(i,k-2,j,2)=f(i,k-2,j,2)-rc*f(i,k,j,2)
      ! cancel part of row d
      rd=d(i,k,j)/b(i,k-2,j)
      d(i,k-2,j)=d(i,k-2,j)-rd*a(i,k-2,j)
      f(i,1,j,1)=f(i,1,j,1)-rd*f(i,k,j,1)
      f(i,1,j,2)=f(i,1,j,2)-rd*f(i,k,j,2)
      ! cancel part of row e
      re=e(i,k,j)/b(i,k-2,j)
      e(i,k-2,j)=e(i,k-2,j)-re*a(i,k-2,j)
      f(i,2,j,1)=f(i,2,j,1)-re*f(i,k,j,1)
      f(i,2,j,2)=f(i,2,j,2)-re*f(i,k,j,2)
    enddo
! same as before for arrays d,e,f, but here no more upper diagonal
    do k=4,3,-1
      ! cancel part of row d
      rd=d(i,k,j)/b(i,k-2,j)
      d(i,k-2,j)=d(i,k-2,j)-rd*a(i,k-2,j)
      f(i,1,j,1)=f(i,1,j,1)-rd*f(i,k,j,1)
      f(i,1,j,2)=f(i,1,j,2)-rd*f(i,k,j,2)
      ! cancel part of row e
      re=e(i,k,j)/b(i,k-2,j)
      e(i,k-2,j)=e(i,k-2,j)-re*a(i,k-2,j)
      f(i,2,j,1)=f(i,2,j,1)-re*f(i,k,j,1)
      f(i,2,j,2)=f(i,2,j,2)-re*f(i,k,j,2)
    enddo
    ! remains only 2x2 array, diagonal and lower diagonal
    !  |d(1), d(2)|  |f(1)|
    !  |e(1), e(2)|  |f(2)|
    ! variable stored in f(1)
    f(i,1,j,1)=(f(i,1,j,1)-f(i,2,j,1)*d(i,2,j)/e(i,2,j))/(d(i,1,j)-e(i,1,j)*d(i,2,j)/e(i,2,j))
    f(i,1,j,2)=(f(i,1,j,2)-f(i,2,j,2)*d(i,2,j)/e(i,2,j))/(d(i,1,j)-e(i,1,j)*d(i,2,j)/e(i,2,j))
    ! variable stored in f(2)
    f(i,2,j,1)=(f(i,2,j,1)-f(i,1,j,1)*e(i,1,j))/(e(i,2,j))
    f(i,2,j,2)=(f(i,2,j,2)-f(i,1,j,2)*e(i,1,j))/(e(i,2,j))

    ! forward-substitute variable from k=2,Nz
    ! solution stored in f
    do k=3,npsiz
      f(i,k,j,1)=(f(i,k,j,1)-a(i,k-2,j)*f(i,k-2,j,1))/b(i,k-2,j)
      f(i,k,j,2)=(f(i,k,j,2)-a(i,k-2,j)*f(i,k-2,j,2))/b(i,k-2,j)
    enddo

  enddo
enddo

return
end
