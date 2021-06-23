      SUBROUTINE polycoeff(f,x,degree,coeff)

!      include 'precision.h'
!      include '../FLW/param.h'
!      include '../DIFFINT/dicommon.h'
      integer degree,i,j
      REAL(8) :: Mc(degree+1,degree+1) 
      REAL(8) :: f(degree+1) 
      REAL(8) :: x(degree+1) 
      REAL(8) :: coeff(degree+1)

      Mc(:,:)=0d0
      coeff(:)=0d0

      do i = 1,int(degree+1),1
      do j = 1,int(degree+1),1
            Mc(i,j)=x(i)**(int(degree+1-j))
!            if (rank .eq. 0) then
!                  print*,'i,j,Mc(i,j)',i,j,Mc(i,j)
!            endif
      enddo
      enddo



      call Gauss(Mc,int(degree+1))

      do i = 1,degree+1
      do j = 1,degree+1
            coeff(i)=coeff(i)+Mc(i,j)*f(j)
      enddo
      enddo

      return
      END






      ! -------------------------------------------------------------------- 
      SUBROUTINE Gauss (a,n)       ! Invert matrix by Gauss method 
      ! -------------------------------------------------------------------- 
      IMPLICIT NONE 
      INTEGER :: n 
      REAL(8) :: a(n,n) 
      ! - - - Local Variables - - - 
      REAL(8) :: b(n,n), c, d, temp(n) 
      INTEGER :: i, j, k, m, imax(1), ipvt(n) 
      ! - - - - - - - - - - - - - - 
      b = a 
      ipvt = (/ (i, i = 1, n) /) 
      DO k = 1,n 
         imax = MAXLOC(ABS(b(k:n,k))) 
         m = k-1+imax(1) 
         IF (m /= k) THEN 
            ipvt( (/m,k/) ) = ipvt( (/k,m/) ) 
            b((/m,k/),:) = b((/k,m/),:) 
         END IF 
         d = 1/b(k,k) 
         temp = b(:,k) 
         DO j = 1, n 
            c = b(k,j)*d 
            b(:,j) = b(:,j)-temp*c 
            b(k,j) = c 
         END DO 
         b(:,k) = temp*(-d) 
         b(k,k) = d 
      END DO 
      a(:,ipvt) = b 
      END SUBROUTINE Gauss 
