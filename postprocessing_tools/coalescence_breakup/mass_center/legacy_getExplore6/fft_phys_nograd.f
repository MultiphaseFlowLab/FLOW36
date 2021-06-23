      subroutine FFT_PHYS_nograd(A,IM,JM,KM,IMH,BU,PHI) 
      
      IMPLICIT NONE
      DOUBLE PRECISION A(IMH,KM,JM,8)
      INTEGER IM,JM,KM,IMH,I,J,K
      DOUBLE PRECISION BU(IM,KM,JM,4), PHI(IM,KM,JM)

      CALL FFTBWDXXINI(A(1,1,1,1),A(1,1,1,1),IM,JM,KM)
      CALL FFTBWDXXINI(A(1,1,1,3),A(1,1,1,3),IM,JM,KM)
      CALL FFTBWDXXINI(A(1,1,1,5),A(1,1,1,5),IM,JM,KM)
      CALL FFTBWDXXINI(A(1,1,1,7),A(1,1,1,7),IM,JM,KM)

C
c
      CALL DCTXX( A(1,1,1,1),A(1,1,1,1), IMH, KM, 2*JM, 2 )
      CALL DCTXX( A(1,1,1,3),A(1,1,1,3), IMH, KM, 2*JM, 2 )
      CALL DCTXX( A(1,1,1,5),A(1,1,1,5), IMH, KM, 2*JM, 2 )
      CALL DCTXX( A(1,1,1,7),A(1,1,1,7), IMH, KM, 2*JM, 2 )  


C
C
C
C     Rearrange arrays 
C     see new version of FFTBWDXX to see that coming from
C     Fourier space in the first half of the arrays there
C     is the first half of the domain, while in Fourier
C     space there was the real part of the quantity.
C
c

! to renormalize the gradient to get the normal field: grad(phi)/|grad(phi)|
      DO 110 J=1,JM
      DO 110 K=1,KM
      DO 110 I=1,IMH-1
      BU(I,K,J,1)         = A(I,K,J,1)
      BU(I+IMH-1,K,J,1)   = A(I,K,J,2)
      BU(I,K,J,2)         = A(I,K,J,3)
      BU(I+IMH-1,K,J,2)   = A(I,K,J,4)
      BU(I,K,J,3)         = A(I,K,J,5)
      BU(I+IMH-1,K,J,3)   = A(I,K,J,6)
      PHI(I,K,J)         = A(I,K,J,7)
      PHI(I+IMH-1,K,J)   = A(I,K,J,8)
  110 CONTINUE      
C
C

      END
