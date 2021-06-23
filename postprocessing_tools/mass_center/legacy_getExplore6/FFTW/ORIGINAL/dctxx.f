      SUBROUTINE DCTXX(X,Y,W,W2,IMH,KM,LM,ITN)
C...Translated by FPP 5.0 (3.03N6) 04/14/92  14:29:20   -dc
C
C **********************************************************************
C *                                                                    *
C *   THIS ROUTINE PERFORMS A DISCRETE CHEBYSHEV TRANSFORM OF LENGTH   *
C *   KM ON IMH*LM VECTORS                                             *
C *                                                                    *
C *   X      -    INPUT VECTOR DIM X(IMH,KM,LM)                        *
C *   Y      -    OUTPUT VECTOR DIM Y(IMH,KM,LM)                       *
C *   W      -    WORK VECTOR DIM W(IMH,KM,LM)                         *
C *   W2     -    WORK VECTOR DIM W(IMH,LM)                            *
C *                                                                    *
C *   ITN    -    1 CHEB. TRANSFORM ( Z TO N )                         *
C *               2 INVERSE CHEB. TRANSFORM ( N TO Z )                 *
C *                                                                    *
C **********************************************************************
C
C     Modified to use FFTW version 3.x
C
      include 'precision.h'
      DIMENSION X(IMH,KM,LM),Y(IMH,KM,LM),W2(1)
      DIMENSION W(KM*2,LM)
      include 'param.h'
      include 'fftw3.f'

      integer*8 plan

      KLEN=KM*2-2

      call dfftw_plan_many_dft_r2c(plan,1,klen,lm,w,klen,1,2*km,
     &           w,klen,1,km,FFTW_MEASURE)

      IF (ITN.EQ.1) then
C
C     Begin forward transform
C
         DO I=1,IMH

            DO K=1,KM
               DO L=1,LM
                  W(K,L)=X(I,K,L)
                  W(K+KM,L)=0.0
               enddo
            enddo

            DO L=1,LM
               W(1,L)=0.5*W(1,L)
               W(KM,L)=0.5*W(KM,L)
            enddo

c            KLEN=KM*2-2

c            call dfftw_plan_many_dft_r2c(plan,1,klen,lm,w,klen,1,2*km,
c     &           w,klen,1,km,FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c(plan,w,w)
      
            DO L=1,LM
               DO K=1,KM
                  Y(I,K,L)=W(2*K-1,L)*2/(KM-1)
               enddo
            enddo

            DO L=1,LM
               Y(I,KM,L)=0.5*Y(I,KM,L)
            enddo
            
c            call dfftw_destroy_plan(plan)

         enddo
         
      else
C
C     Begin Inverse Transform
C
         DO I=1,IMH

            DO K=1,KM
               DO L=1,LM
                  W(K,L)=X(I,K,L)
                  W(K+KM,L)=0.0
               enddo
            enddo

            DO L=1,LM
               W(1,L)=0.5*W(1,L)
            enddo

c            KLEN=2*KM-2

c            call dfftw_plan_many_dft_r2c(plan,1,klen,lm,w,klen,1,2*km,
c     &           w,klen,1,km,FFTW_ESTIMATE)
            call dfftw_execute_dft_r2c(plan,w,w)

            DO L=1,LM
               DO K=1,KM
                  Y(I,K,L)=W(2*K-1,L)
               enddo
            enddo

c            call dfftw_destroy_plan(plan)

         enddo

      endif

      call dfftw_destroy_plan(plan) 
    
      RETURN
      END
