      SUBROUTINE ACGEN(CA,KSET,CFILE)

c-----------------------------------------------------------------------

      IMPLICIT NON E 
      INTEGER KSET, ICA 
      CHARACTER*5 CFILE
      REAL*8 CA(15)

c-----------------------------------------------------------------------


      IF ((CA(KSET).GE.1.0).AND.(CA(KSET).LE.10)) THEN
          ICA=CA(KSET)*1000.0
          WRITE(CFILE(2:5),'(i4)') ICA
          CFILE(1:1)=CFILE(2:2)
          CFILE(2:2)='_'
      ENDIF

      IF ((CA(KSET).GE.0.1).AND.(CA(KSET).LE.1.0)) THEN
          ICA=CA(KSET)*1000.0
          WRITE(CFILE(3:5),'(i3)') ICA
          CFILE(1:1)='0'
          CFILE(2:2)='_'
      ENDIF     


c-----------------------------------------------------------------------
      RETURN
      END
