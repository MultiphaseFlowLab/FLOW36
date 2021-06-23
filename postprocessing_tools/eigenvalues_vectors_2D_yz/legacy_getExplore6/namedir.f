      SUBROUTINE NAMEDIR(STRNDIR,IDIR,NLDIRMAX)
      
      integer IDIR,NLDIRMAX,i
      character (len=NLDIRMAX) STRNDIR

      write(*,*)'dentro namedi'

      IDIR=0
   10 IDIR=IDIR+1
      if(STRNDIR(IDIR:IDIR).eq.' ') goto 10
      if(IDIR.gt.1) then
       do i=1,NLDIRMAX+1-IDIR
        STRNDIR(i:i)=STRNDIR(i+IDIR-1:i+IDIR-1)
       enddo
       do i=NLDIRMAX-IDIR,NLDIRMAX
        STRNDIR(i:i)=' '
       enddo
      end if

      IDIR=NLDIRMAX+1
   11 IDIR=IDIR-1
      if(STRNDIR(IDIR:IDIR).eq.' ') goto 11

      if(STRNDIR(IDIR:IDIR).eq.'/') then
       IDIR=IDIR-1
      end if

      do I=1,IDIR
       if(STRNDIR(I:I).eq.' ') then
        print *, 'ERRORE in file di input !!!'
        print *, 'Nome della directory =', STRNDIR(1:IDIR)
        print *, ' NON deve contenere spazi !!!'
        STOP
       end if
      enddo

      write(*,*)'in fondo namedir'

      RETURN
      END

