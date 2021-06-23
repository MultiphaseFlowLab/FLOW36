      SUBROUTINE READFT12_1(KSET,DIRINP,IDIRINP,FINAME1,IMH,JMP,
     &                  NX,NY,NZ,NNT,A,NTIM,TIME) 

c-----------------------------------------------------------------------

      IMPLICIT NONE 
      INTEGER KSET,KNODE
      INTEGER I1,J1,J2,I,J,K
      INTEGER IDIRINP
      INTEGER IMH,JMP,IFLD,NX,NY,NZ
      CHARACTER*26 FINAME1
      CHARACTER*85 DIRINP
      CHARACTER*10 TFILE
      INTEGER NNT, UNT

      DOUBLE PRECISION A(IMH,NZ,NY,8)
      DOUBLE PRECISION TIME
      INTEGER NTIM, ITFILE

c-----------------------------------------------------------------------

      

      IF ((KSET.GE.0).AND.(KSET.LE.9)) THEN
          WRITE(FINAME1(26:26),'(i1)') KSET
          FINAME1(24:25)='00'
      ENDIF

      IF ((KSET.GE.10).AND.(KSET.LE.99)) THEN
          WRITE(FINAME1(25:26),'(i2)') KSET
          FINAME1(24:24)='0'
      ENDIF

      IF ((KSET.GE.100).AND.(KSET.LE.999)) THEN
          WRITE(FINAME1(24:26),'(i3)') KSET
      ENDIF

      UNT=KSET+101

c      OPEN(unit=UNT,FILE=DIRINP(1:IDIRINP)//'/'//TFILE(1:ITFILE)//
c     $     '/OUTPUT/'//FINAME1,STATUS='OLD',FORM='UNFORMATTED')

c read ft12_diffintXXXXXX-rankXXX if NNT > 1, otherwis
c read ft12_diffintXXXXXX
      IF (NNT.GT.1) THEN

         OPEN(unit=UNT,FILE=DIRINP(1:IDIRINP)//'/'//TFILE(1:ITFILE)//
     $        FINAME1,STATUS='OLD',FORM='UNFORMATTED')     

         WRITE(*,*)'Reading:',DIRINP(1:IDIRINP)//'/'//TFILE(1:ITFILE)//
     $        FINAME1
      ELSE

         OPEN(unit=UNT,FILE=DIRINP(1:IDIRINP)//'/'//TFILE(1:ITFILE)//
     $        FINAME1(1:18),STATUS='OLD',FORM='UNFORMATTED')

         WRITE(*,*)'Reading:',DIRINP(1:IDIRINP)//'/'//TFILE(1:ITFILE)//
     $        FINAME1(1:18)
      ENDIF
     
C - reading of the file ft12_diffint...

      READ (UNT) NTIM, TIME

      I1=NX/3+1
      J1=NY/3+1
      J2=NY+2-J1
      DO 10 IFLD=1,8
      DO 10 I=1,I1
      DO 10 J=1,JMP
c
      IF ((J+KSET*JMP).GT.J1 .AND. (J+KSET*JMP).LT.J2) GO TO 10
            READ(UNT) (A(I,K,(J+KSET*JMP),IFLD), K=1,NZ)
 10   CONTINUE

  

c-----------------------------------------------------------------------
      RETURN
      END
