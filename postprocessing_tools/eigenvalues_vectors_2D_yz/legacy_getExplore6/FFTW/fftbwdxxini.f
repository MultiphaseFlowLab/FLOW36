      SUBROUTINE FFTBWDXXINI(X,Y,IM,JM,KM)
      include 'precision.h'
      include '../param.h'
      include 'fftw3.f'
C
C     Modified to use FFTW version 3.x
C
C - INI version of FFTBWDXX collects the physical space variables  
C - on output array Y considering Y(:,:,:,1) and Y(:,:,:,2) as the 
C - first half and the second half of the domain in the length-wise 
C - direction.
C - Only the last do-loop diiffers from FFTBWDXX version.
C - L. Scarbolo May 2012
C  
      DIMENSION X(IM/2+1,KM,JM,2),Y(IM/2+1,KM,JM,2)

C - L.Scarbolo May 2012
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: W
      integer*8 plan

      allocate(W(IM+2,KM,JM))

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
C - do-loop variations
      DO J=1,JM
         DO K=1,KM
            DO I=1,IMH-1
               Y(I,K,J,1)=W(I,K,J)/(IM*JM)
               Y(I,K,J,2)=W(IMH+I-1,K,J)/(IM*JM)
            enddo
         enddo
      enddo
     
      deallocate(W)
 
      call dfftw_destroy_plan(plan)

      RETURN
      END
