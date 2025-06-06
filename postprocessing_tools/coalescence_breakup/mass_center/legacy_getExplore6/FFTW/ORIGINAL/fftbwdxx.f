      SUBROUTINE FFTBWDXX(X,Y,W,IM,JM,KM)
      include 'precision.h'
      include 'param.h'
      include 'fftw3.f'
C
C     Modified to use FFTW version 3.x
C
      DIMENSION X(IM/2+1,KM,JM,2),Y(IM/2+1,KM,JM,2),W(IM+2,KM,JM)

      integer*8 plan

      IMH=IM/2+1
      
      DO J=1,JM
         DO K=1,KM
            DO I=1,IMH
               W(I*2-1,K,J)=X(I,K,J,1)
               W(I*2,K,J)=X(I,K,J,2)
            enddo
         enddo
      enddo

      call dfftw_plan_many_dft(plan,1,jm,imh*km,w,jm,imh*km,1,w,jm,
     &     imh*km,1,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,w,w)

      call dfftw_plan_many_dft_c2r(plan,1,im,jm*km,w,im,1,imh,w,im,
     &     1,im+2,FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(plan,w,w)

      DO J=1,JM
         DO K=1,KM
            DO I=1,IMH-1
               Y(I,K,J,1)=W(I*2-1,K,J)/(IM*JM)
               Y(I,K,J,2)=W(I*2,K,J)/(IM*JM)
            enddo
         enddo
      enddo
 
      call dfftw_destroy_plan(plan)

      RETURN
      END
