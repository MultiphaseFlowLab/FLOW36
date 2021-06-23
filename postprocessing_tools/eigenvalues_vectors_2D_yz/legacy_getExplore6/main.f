c***********************************************************************
c  FORPLOT
C  Program to create file for plotting.
c  Different formats are included 
C c***********************************************************************
      PROGRAM MAIN 

C -1- C VARIABLES DECLARATION

      IMPLICIT NONE

      include 'param.h'
      
      INTEGER KNODE,NODENO
      CHARACTER*26 FINAME1     
      CHARACTER*19 FINAME2
      CHARACTER*16 FINAME3
      INTEGER IT, STEPIN, STEPFIN, DSTEP
      INTEGER K, I, J, ip 
      CHARACTER*85 DIRINP, DIROUT,PDIRINP
      CHARACTER*6 AFILE, AFILE2
      INTEGER IDIRINP, IDIROUT, ITFILE,IPDIRINP
      REAL*8 LX, LY, LZ, RE, CA, CH 
      REAL*8 X(NX), Y(NY), Z(NZ)   
      REAL*8 DX, DY, DZ, PI, CC
      DOUBLE PRECISION TIME
      INTEGER NTIM
      INTEGER KK, NX1, NX2, NY1, NY2, RANK, DNZ,DNX,DNY
      DOUBLE PRECISION interface_coord(6,nx*ny*nz)
      INTEGER NP_I,things_to_do(4)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: A
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHI


C -2- READING OF INPUT FILE
    
      CALL GETINPUTDATA(STEPIN, STEPFIN, DSTEP,DIRINP,DIROUT, LX, LY, LZ, RE, CA, CH, NX1,NY1,NX2,NY2,DNZ,DNX,DNY,things_to_do)
      

      WRITE(*,*) 'AFTER GETINPUTDATA'

      CALL NAMEDIR(DIROUT,IDIROUT,len(DIROUT))
      CALL NAMEDIR(DIRINP,IDIRINP,len(DIRINP)) 


      WRITE(*,*) 'stepin:',STEPIN

C - calculation of constant quatities

      PI = 4.0*ATAN(1.0)
      PRINT*,PI

      DO K = 1, NZ
          Z(K)=COS((DBLE(K-1)*PI)/DBLE(NZ-1))
      ENDDO

      DX = (LX/(NX-1))
      DY = (LY/(NY-1))

      X(1) = 0.0
      Y(1) = 0.0

      DO I = 2, NX
         X(I) = X(I-1) + DX
      ENDDO

      DO J = 2, NY
         Y(J) = Y(J-1) + DY
      ENDDO
      
      
      FINAME1(1:12)='ft12_diffint'
      FINAME1(19:24)='-rank'
      FINAME2(1:9)='OUT_FEDBI'
      FINAME2(16:19)='.vtk'
c      FINAME3(1:10)='ft12_press'
      OPEN(11,FILE=DIROUT(1:IDIROUT)//'/Drops.txt',POSITION='APPEND',STATUS='NEW')
      write(11,*) '% This file contains:'
      write(11,*) '% Time,Drops number, Mean Volume, Equivalent Radius, Mean Volume(Small drops),Reqm(Small Drops)'            
      close(11) 
      
      
C - Time step advancement
      CC=0.0
      
      DO 201 IT=STEPIN,STEPFIN,DSTEP
        
         CC = CC + 1.0

         WRITE(*,*) ' '
         WRITE(*,*) '************    TIME STEP:' ,IT, '****************'
         WRITE(*,*) 'Starting to analyze'
 
         CALL AFGEN(IT,AFILE,AFILE2)

         FINAME1(13:18) = AFILE
         FINAME2(10:15) = AFILE2
         FINAME3(11:16) = AFILE

c - Read input files - FFT coefficients
         ALLOCATE(A(NXH,NZ,NY,8))
         WRITE(*,*) 'Setting A(:,:,:,:) to ZERO!!!'
         A(:,:,:,:)=0d0
         DO RANK=0,(NNT-1)
         WRITE(*,*)'--------- READ ft12, RANK:',RANK,'-----------'
              CALL READFT12_1(RANK,DIRINP,IDIRINP,FINAME1,NXH,NYP,NX,NY,NZ,NNT,A,NTIM,TIME)
         ENDDO

         WRITE(*,*) 'File ',FINAME1(1:18),'has been red, now fftbwd'  
         ALLOCATE(PHI(NX,NZ,NY))
         PHI(:,:,:)=0.0d0
         print*,'Call fftw'
         CALL FFT_PHYS(A,NX,NY,NZ,NXH,PHI)
         DEALLOCATE(A)
         print*,'Call Explore'
         CALL EXPLORE(PHI,NX,NY,NZ,NXH,X,Y,Z,DX,DY,TIME,DIROUT,IDIROUT)
         
        DEALLOCATE(PHI)         

 201  CONTINUE 
 
 1001 FORMAT(A26)
 1002 FORMAT(A33)
 1003 FORMAT(A5)
 1004 FORMAT(A24)
 1005 FORMAT(A10,3(1X,I3))
 1006 FORMAT(A13,1X,I3,1X,A6)
 1007 FORMAT(E12.6)
 1008 FORMAT(A10,1X,I7)
 1009 FORMAT(A7,1X,A5,1X,A6,1X,I1)
 1010 FORMAT(A20)
 1011 FORMAT(E16.8)
 
      write(*,*)
      write(*,*)'********************************************'
      write(*,*)'--------- THIS IS THE END ------------------'
      write(*,*)'********************************************'



      END

c======================================================================
