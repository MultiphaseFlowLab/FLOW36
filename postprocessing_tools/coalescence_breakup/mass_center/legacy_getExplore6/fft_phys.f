      subroutine FFT_PHYS(A,IM,JM,KM,IMH,PHI) 
      
      IMPLICIT NONE
      DOUBLE PRECISION A(IMH,KM,JM,8)
      INTEGER IM,JM,KM,IMH,I,J,K
      DOUBLE PRECISION PHI(IM,KM,JM)

      CALL FFTBWDXXINI(A(1,1,1,7),A(1,1,1,7),IM,JM,KM)

      CALL DCTXX( A(1,1,1,7),A(1,1,1,7), IMH, KM, 2*JM, 2 )  


      DO 110 J=1,JM
      DO 110 K=1,KM
      DO 110 I=1,IMH-1
      PHI(I,K,J)          = A(I,K,J,7)
      PHI(I+IMH-1,K,J)    = A(I,K,J,8)
  110 CONTINUE      
      END
