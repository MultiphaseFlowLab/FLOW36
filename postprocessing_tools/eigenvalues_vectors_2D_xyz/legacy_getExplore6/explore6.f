      SUBROUTINE EXPLORE(PHI,NX,NY,NZ,NXH,X,Y,Z,DX,DY,TIME,DIROUT,IDIROUT)
c-----------------------------------------------------------------------
C-----A.Roccon 07/2015
C-----Subroutine explore the PHI field....and then counts the drop and calculate the medimum mass of the drops

      IMPLICIT NONE 
      INTEGER NX, NXH, NY, NZ, I, J, K, II, KK, JJ, L, S, T, F
      INTEGER DROPS,SDROPS,IM,IP,JM,JP
      DOUBLE PRECISION X(NX), Y(NY), Z(NZ) 
      DOUBLE PRECISION DV,VOL,VOLM,REQ,REQ2,REQM,LIM,VOLM2,REQM2
      DOUBLE PRECISION PHI (NX,NZ,NY)
      INTEGER IDIROUT
      real*8 time
      DOUBLE PRECISION DX, DY, DZ, PI
      CHARACTER*85 DIROUT
      CHARACTER*6 AFILE2
      CHARACTER*132 filename
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHIE
C      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MAP
      INTEGER P(80000,3), PN(80000,3)
      
      PI = 4.0*ATAN(1.0)
      
         
      ALLOCATE(PHIE(NX,NZ,NY))

      
      DO I=1,NX
      DO J=1,NY
      DO K=1,NZ
                  PHIE(I,K,J)=PHI(I,K,J)
      ENDDO
      ENDDO
      ENDDO
      
C      print*,'PHIE copied....'
      print*,'Start the Dora-Exploration loop'
      print*,'New version'
      DROPS=0
      VOLM=0.0d0
      VOLM2=0.0d0
      REQM=0.0d0
      REQM2=0.0d0
      P(:,:)=0
      DO KK=2,NZ-1
      DO JJ=1,NY
      DO II=1,NX
          IF (PHIE(II,KK,JJ) .GT. 0.0d0) THEN
              PHIE(II,KK,JJ)=-1.0d0
              DROPS=DROPS+1
              DZ=0.5d0*ABS((Z(KK-1)-Z(KK+1)))
              VOL=DX*DY*DZ
              P(1,1)=II
              P(1,2)=KK
              P(1,3)=JJ
              DO S=1,400
              F=1
              PN(:,:)=0
              DO L=1,80000 !(it's enough???)
                  IF( P(L,1).GT.0) THEN !! check if the P-list is valid or not
                  I=P(L,1)
                  k=P(L,2)
                  J=P(L,3)
                      DZ=0.5d0*ABS((Z(K-1)-Z(K+1)))
                          VOL=VOL+DX*DY*DZ
                          IP=I+1 ! I plus 1
                          IM=I-1
                          JP=J+1
                          JM=J-1
                          IF (IP.EQ.NX+1) THEN
                              IP=1
c                              print*,'ip'
                          ENDIF
                          IF (IM.EQ.0) THEN
                              IM=NX
C                              print*,'im'
                          ENDIF
                          IF (JP.EQ.NY+1) THEN
                              JP=1
C                              print*,'jp'
                          ENDIF
                          IF (JM.EQ.0) THEN
                              JM=NY
C                              print*,'jm'
                          ENDIF
C                          print*,IM,K-1,J
                          IF (PHIE(IM,K-1,J).GT.0.0d0)  THEN
                              PHIE(IM,K-1,J)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K-1
                              PN(F,3)=J
                              F=F+1
C                             VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IM,K-1,JP).GT.0.0d0)  THEN
                              PHIE(IM,K-1,JP)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K-1
                              PN(F,3)=JP
                              F=F+1
                          ENDIF
                          IF (PHIE(IM,K-1,JM).GT.0.0d0 ) THEN
                              PHIE(IM,K-1,JM)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K-1
                              PN(F,3)=JM
                              F=F+1
                          ENDIF
                          !!!!!!! PLANE I-1-----K
                          IF (PHIE(IM,K,J).GT.0.0d0 )  THEN
                              PHIE(IM,K,J)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K
                              PN(F,3)=J
                              F=F+1
                          ENDIF
                          IF (PHIE(IM,K,JP).GT.0.0d0 )  THEN
                              PHIE(IM,K,JP)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K
                              PN(F,3)=JP
                              F=F+1
                              
                          ENDIF
                          IF (PHIE(IM,K,JM).GT.0.0d0 )  THEN
                              PHIE(IM,K,JM)=-1.0d0 
                            !  MAP(IM,K,JM)=1
                              PN(F,1)=IM
                              PN(F,2)=K
                              PN(F,3)=JM
                              F=F+1
                          !    VOL=VOL+DX*DY*DZ
                              
                          ENDIF
                          !!!!!!! PLANE IM-----K+1
                          IF (PHIE(IM,K+1,J).GT.0.0d0 )  THEN
                              PHIE(IM,K+1,J)=-1.0d0 
                           !   MAP(IM,K+1,J)=1
                              PN(F,1)=IM
                              PN(F,2)=K+1
                              PN(F,3)=J
                              F=F+1
                          !    VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IM,K+1,JP).GT.0.0d0 )  THEN
                              PHIE(IM,K+1,JP)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K+1
                              PN(F,3)=JP
                              F=F+1
                           !   VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IM,K+1,JM).GT.0.0d0 )  THEN
                              PHIE(IM,K+1,JM)=-1.0d0 
                              PN(F,1)=IM
                              PN(F,2)=K+1
                              PN(F,3)=JM
                              F=F+1
                           !   VOL=VOL+DX*DY*DZ
                          ENDIF
                          !!!!PLANE I
                          !!! PLANE I K_1
                          IF (PHIE(I,K-1,J).GT.0.0d0 )  THEN
                              PHIE(I,K-1,J)=-1.0d0 
                              PN(F,1)=I
                              PN(F,2)=K-1
                              PN(F,3)=J
                              F=F+1
                           !   VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(I,K-1,JP).GT.0.0d0 )  THEN
                              PHIE(I,K-1,JP)=-1.0d0 
                           !   MAP(I,K-1,JP)=1
                              PN(F,1)=I
                              PN(F,2)=K-1
                              PN(F,3)=JP
                              F=F+1
                           !   VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(I,K-1,JM).GT.0.0d0 )  THEN
                              PHIE(I,K-1,JM)=-1.0d0 
                              PN(F,1)=I
                              PN(F,2)=K-1
                              PN(F,3)=JM
                              F=F+1
                         !     VOL=VOL+DX*DY*DZ
                              
                          ENDIF
                          !!!!!!! PLANE I-----K
                          IF (PHIE(I,K,JP).GT.0.0d0 )  THEN
                              PHIE(I,K,JP)=-1.0d0 
                              PN(F,1)=I
                              PN(F,2)=K
                              PN(F,3)=JP
                              F=F+1
                         !     VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(I,K,JM).GT.0.0d0 )  THEN
                              PHIE(I,K,JM)=-1.0d0 
                           !   MAP(I,K,JM)=1
                              PN(F,1)=I
                              PN(F,2)=K
                              PN(F,3)=JM
                              F=F+1
                          !    VOL=VOL+DX*DY*DZ
                          ENDIF
                          !!!!!!! PLANE IM-----K+1
                          IF (PHIE(I,K+1,J).GT.0.0d0 )  THEN
                              PHIE(I,K+1,J)=-1.0d0 
                           !   MAP(I,K+1,J)=1
                              PN(F,1)=I
                              PN(F,2)=K+1
                              PN(F,3)=J
                              F=F+1
                         !     VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(I,K+1,JP).GT.0.0d0 )  THEN
                              PHIE(I,K+1,JP)=-1.0d0 
                        !      MAP(I,K+1,JP)=1
                              PN(F,1)=I
                              PN(F,2)=K+1
                              PN(F,3)=JP
                              F=F+1
                          !    VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(I,K+1,JM).GT.0.0d0 )  THEN
                              PHIE(I,K+1,JM)=-1.0d0 
                         !     MAP(I,K+1,JM)=1
                              PN(F,1)=I
                              PN(F,2)=K+1
                              PN(F,3)=JM
                              F=F+1
                          !    VOL=VOL+DX*DY*DZ
                          ENDIF
                          !!!PLANE IP
                          !!!!PLANE IP---K-1
                          IF (PHIE(IP,K-1,J).GT.0.0d0 )  THEN
                              PHIE(IP,K-1,J)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K-1
                              PN(F,3)=J
                              F=F+1
                         !     VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IP,K-1,JP).GT.0.0d0 )  THEN
                              PHIE(IP,K-1,JP)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K-1
                              PN(F,3)=JP
                              F=F+1
                       !       VOL=VOL+DX*DY*DZ
                          ENDIF                          
                          IF (PHIE(IP,K-1,JM).GT.0.0d0 )  THEN
                              PHIE(IP,K-1,JM)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K-1
                              PN(F,3)=JM
                              F=F+1
                        !      VOL=VOL+DX*DY*DZ
                          ENDIF
                          !!!!!!! PLANE IP-----K
                          IF (PHIE(IP,K,J).GT.0.0d0 )  THEN
                              PHIE(IP,K,J)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K
                              PN(F,3)=J
                              F=F+1
                         !     VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IP,K,JP).GT.0.0d0 )  THEN
                              PHIE(IP,K,JP)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K
                              PN(F,3)=JP
                              F=F+1
                          !    VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IP,K,JM).GT.0.0d0 )  THEN
                              PHIE(IP,K,JM)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K
                              PN(F,3)=JM
                              F=F+1
                           !   VOL=VOL+DX*DY*DZ
                          ENDIF
                          !!!!!!! PLANE IP-----K+1
                          IF (PHIE(IP,K+1,J).GT.0.0d0 )  THEN
                              PHIE(IP,K+1,J)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K+1
                              PN(F,3)=J
                              F=F+1
                            !  VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IP,K+1,JP).GT.0.0d0 )  THEN
                              PHIE(IP,K+1,JP)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K+1
                              PN(F,3)=JP
                              F=F+1
                             ! VOL=VOL+DX*DY*DZ
                          ENDIF
                          IF (PHIE(IP,K+1,JM).GT.0.0d0 )  THEN
                              PHIE(IP,K+1,JM)=-1.0d0 
                              PN(F,1)=IP
                              PN(F,2)=K+1
                              PN(F,3)=JM
                              F=F+1
                        !      VOL=VOL+DX*DY*DZ
                          ENDIF                          
                  ENDIF
              ENDDO
C              print*,'travaso PN-->P'
              DO T=1,80000
                  P(T,1)=PN(T,1)
                  P(T,2)=PN(T,2)
                  P(T,3)=PN(T,3)
              ENDDO
C              print*,'erase PN'
              PN(:,:)=0
              ENDDO
C              print*,'exiting the loop of exploration of the n-th drop'
              LIM=4.0d0
              IF (VOL.LT.LIM) THEN
                  REQ2=(4*VOL/(3*PI))**(1.0d0/3.0d0) !! RAGGIO EQUIVALENTE GOCCIA SE NON é TROPPO GRANDE...
                  SDROPS=SDROPS+1
                  VOLM2 = (VOLM2*(SDROPS-1)+VOL)/(SDROPS)
                  REQM2 = (REQM2*(SDROPS-1)+REQ2)/(SDROPS)
              ENDIF
              REQ=(4*VOL/(3*PI))**(1.0d0/3.0d0)
C              print*,'DROP N',drops,'REQ',REQ,'VOL',VOL
              VOLM = (VOLM*(DROPS-1)+VOL)/(DROPS)
              REQM = (REQM*(DROPS-1)+REQ)/(DROPS)
              VOL=0.0d0
              REQ=0.0d0
              REQ2=0.0d0
          ENDIF
        ENDDO
        ENDDO
        ENDDO
        ! FINE DI TUTTA L'ESPLOTAZIONE...
        !REQ=(4*VOLM/(3*PI))**(1.0d0/3.0d0)
        print*,'volm',VOLM
        print*,'Rad.eq',REQM
        print*,'drops',DROPS
        print*,'volm',VOLM2
        print*,'Rad.eq',REQM2
        print*,'drops',SDROPS
      
                  
      
         OPEN(56,FILE=DIROUT(1:IDIROUT)//'/Drops.txt',POSITION='APPEND')
         write(56,1007) TIME,DROPS,SDROPS,VOLM,REQM,VOLM2,REQM2
         close(56)
 1007 FORMAT(E14.8,3x,I8,3x,I8,3x,E14.8,3x,E14.8,3x,E14.8,3x,E14.8)
         
         
         DEALLOCATE(PHIE)


c----------------------------------------------------------------------
      RETURN
      END
