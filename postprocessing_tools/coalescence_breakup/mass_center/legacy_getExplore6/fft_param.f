      subroutine FFT_PARAM(PI,AL1,AL2,IM,JM,IMH,AK1,AK2) 
      
      IMPLICIT NONE
      INTEGER IM,JM,IMH,IMP2,I,J,JMH
      DOUBLE PRECISION AK1(IM), AK2(JM)
      DOUBLE PRECISION PI,AL1,AL2,twopi,waveno
      
      twopi = 2.*PI 
      imp2=im+2
      JMH=JM/2+1

      do 10 i=2,imh
      waveno = (i-1)*twopi/al1
      ak1(i) = waveno
      ak1(imp2-i) = -waveno
   10 continue

      ak1(1) = 0.
      ak2(1) = 0.

      do 20 j=2,jmh
      waveno = (j-1)*twopi/al2
      ak2(j) = waveno
      ak2(jm+2-j) = -waveno
   20 continue

      RETURN
      END

