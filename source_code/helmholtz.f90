subroutine helmholtz(f,beta2,p,q,r,z_p)
! solve equation
! y"-beta2*y=f
! with boundary conditions
! p(i)*y(z(i))+q(i)*y'(z(i))=r(i)       @ z(i) , i=1,2

use commondata
use par_size

double precision, dimension(spx,nz,spy,2) :: f,h
double precision, dimension(spx,spy) :: beta2
double precision, dimension(2) :: p,q,r,z_p
double precision :: a(spx,nz-2,spy),b(spx,nz-2,spy),c(spx,nz-4,spy),d(spx,nz,spy),e(spx,nz,spy)
double precision :: t0,dt0,t,dert

integer :: i,j,k

! system has line 1 and 2 full and then it is tridiagonal from line 3 to Nz
! a is the lower diagonal, b the diagonal, c the upper diagonal, d the first row, e the second row
! in the diagonal are excluded the first 2 rows (already in the arrays d and e)


! assemble diagonals
!$acc kernels
do j=1,spy
  do k=3,nz
    do i=1,spx
      a(i,k-2,j)=-beta2(i,j)*dble(k)
      b(i,k-2,j)=4.0d0*dble(k*(k-1)*(k-2))+dble(2*(k-1))*beta2(i,j)
    enddo
  enddo
enddo
!$acc end kernels

!$acc kernels
do j=1,spy
  do k=3,nz-2
    do i=1,spx
      c(i,k-2,j)=-beta2(i,j)*dble(k-2)
    enddo
  enddo
enddo
!$acc end kernels

! assemble 1st row
!$acc kernels
t0=1.0d0
dt0=0.0d0
d(:,1,:)=0.5d0*(p(1)*t0+q(1)*dt0)
do k=2,nz
! see Canuto et al. 2006, pag. 85
  t=dble((z_p(1))**(k-1))
  dert=dble((z_p(1))**(k)*(k-1)**2)
  d(:,k,:)=p(1)*t+q(1)*dert
enddo
!$acc end kernels

! assemble 2nd ro
!$acc kernels
t0=1.0d0
dt0=0.0d0
e(:,1,:)=0.5d0*(p(2)*t0+q(2)*dt0)
do k=2,nz
! see Canuto et al. 2006, pag. 85
  t=dble((z_p(2))**(k-1))
  dert=dble((z_p(2))**(k)*(k-1)**2)
  e(:,k,:)=p(2)*t+q(2)*dert
enddo
!$acc end kernels

! assemble RHS of equation
!$acc kernels
h(:,1,:,:)=r(1)
h(:,2,:,:)=r(2)

do k=3,nz-2
  h(:,k,:,:)=dble(k)*f(:,k-2,:,:)-dble(2*(k-1))*f(:,k,:,:)+dble(k-2)*f(:,k+2,:,:)
enddo
h(:,nz-1,:,:)=dble(nz-1)*f(:,nz-3,:,:)-dble(2*(nz-2))*f(:,nz-1,:,:)
h(:,nz,:,:)=dble(nz)*f(:,nz-2,:,:)-dble(2*(nz-1))*f(:,nz,:,:)
!$acc end kernels

call gauss_solver(a,b,c,d,e,h)

!$acc kernels
f=h
!$acc end kernels

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
do j=1,spypsi
  do k=3,npsiz
    do i=1,spxpsi
      a(i,k-2,j)=-beta2(i,j)*dble(k)
      b(i,k-2,j)=4.0d0*dble(k*(k-1)*(k-2))+dble(2*(k-1))*beta2(i,j)
    enddo
  enddo
enddo

do j=1,spypsi
  do k=3,nz-2
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
