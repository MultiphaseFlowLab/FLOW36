      SUBROUTINE FFTFWDXX(X,Y,IM,JM,KM)
      include 'precision.h'
      include 'param.h'
      include 'fftw3.f'
C
C     Modified to use FFTW version 3.x
C
      DIMENSION X(IM/2+1,KM,JM,2),Y(IM/2+1,KM,JM,2)

C - L.Scarbolo May 2012
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: W

      integer*8 plan

      allocate(W(IM+2,KM,JM))

      IMH=IM/2+1

      DO K=1,KM
         DO J=1,JM
            DO I=1,IMH-1
               W(2*I-1,K,J)=X(I,K,J,1)
               W(2*I,K,J)=X(I,K,J,2)
            enddo
         enddo
      enddo


      call dfftw_plan_many_dft_r2c(plan,1,im,jm*km,w,im,1,im+2,w,im,
     &     1,imh,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,w,w)

      call dfftw_plan_many_dft(plan,1,jm,imh*km,w,jm,imh*km,1,w,
     &     jm,imh*km,1,FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,w,w)

      DO J=1,JM
         DO K=1,KM
            DO I=1,IMH
               Y(I,K,J,1)=W(2*I-1,K,J)
               Y(I,K,J,2)=W(2*I,K,J)
            enddo
         enddo
      enddo

      deallocate(W)

      call dfftw_destroy_plan(plan)

      RETURN
      END
