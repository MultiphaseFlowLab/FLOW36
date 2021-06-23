      SUBROUTINE AFGEN(IT,AFILE,AFILE2)

c-----------------------------------------------------------------------

      IMPLICIT NONE 
      INTEGER IT 
      CHARACTER*6 AFILE, AFILE2

c-----------------------------------------------------------------------

      IF ((IT.GE.0).AND.(IT.LE.9)) THEN
         AFILE(1:5)='-----'
         AFILE2(1:5)='00000'
         WRITE(AFILE(6:6),'(i1)') IT
         WRITE(AFILE2(6:6),'(i1)') IT
      ENDIF

      IF ((IT.GE.10).AND.(IT.LE.99)) THEN
         AFILE(1:4)='----'
         AFILE2(1:4)='0000'
         WRITE(AFILE(5:6),'(i2)') IT
         WRITE(AFILE2(5:6),'(i2)') IT
      ENDIF

      IF ((IT.GE.100).AND.(IT.LE.999)) THEN
         AFILE(1:3)='---'
         AFILE2(1:3)='000'
         WRITE(AFILE(4:6),'(i3)') IT
         WRITE(AFILE2(4:6),'(i3)') IT
      ENDIF

      IF ((IT.GE.1000).AND.(IT.LE.9999)) THEN
         AFILE(1:2)='--'
         AFILE2(1:2)='00'
         WRITE(AFILE(3:6),'(i4)') IT
         WRITE(AFILE2(3:6),'(i4)') IT
      ENDIF

      IF ((IT.GE.10000).AND.(IT.LE.99999)) THEN
         AFILE(1:1)='-'
         AFILE2(1:1)='0'
         WRITE(AFILE(2:6),'(i5)') IT
         WRITE(AFILE2(2:6),'(i5)') IT
      ENDIF

      IF ((IT.GE.100000).AND.(IT.LE.999999)) THEN
         WRITE(AFILE,'(i6)') IT
         WRITE(AFILE2,'(i6)') IT
      ENDIF

c-----------------------------------------------------------------------
      RETURN
      END
