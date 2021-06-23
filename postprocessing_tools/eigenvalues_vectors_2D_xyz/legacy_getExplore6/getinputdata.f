      SUBROUTINE GETINPUTDATA(STEPIN, STEPFIN, DSTEP,
     &              DIRINP, DIROUT,LX, LY, LZ, RE, 
     &              CA, CH,NX1,NY1,NX2,NY2,
     &              DNZ, DNX,DNY,things_to_do)


c-----------------------------------------------------------------------

      IMPLICIT NONE 
      INTEGER STEPIN, STEPFIN, DSTEP
      REAL*8 LX, LY, LZ, RE, CA, CH 
      CHARACTER*85 DIRINP, DIROUT
      INTEGER NX1,NY1,NX2,NY2,DNZ, DNX,DNY,things_to_do(4)      

c-----------------------------------------------------------------------





      OPEN(UNIT=9,FILE='input.inp',STATUS='OLD',FORM='FORMATTED')
      read(9,*) 
      read(9,*)
      read(9,*) STEPIN
      read(9,*) STEPFIN
      read(9,*) DSTEP 
      read(9,*) 
      read(9,*)
      read(9,*) NX1 
      read(9,*) NY1 
      read(9,*) NX2  
      read(9,*) NY2 
      read(9,*) DNZ 
      read(9,*) DNX 
      read(9,*) DNY 
      read(9,*)
      read(9,*)
      read(9,*) LX 
      read(9,*) LY
      read(9,*) LZ
      read(9,*)
      read(9,*)
      read(9,*) RE
      read(9,*) CA
      read(9,*) CH
      read(9,*)
      read(9,*)
      read(9,8085) DIRINP
      read(9,*)
      read(9,*)
      read(9,8085) DIROUT
      read(9,*)
      read(9,*)
      read(9,*) things_to_do(1)
      read(9,*) things_to_do(2)
      read(9,*) things_to_do(3)
      read(9,*) things_to_do(4)
      read(9,*) 
      read(9,*) 
      CLOSE(UNIT=9,STATUS='KEEP')
 8085 format(A85)

c-----------------------------------------------------------------------
      RETURN
      END
