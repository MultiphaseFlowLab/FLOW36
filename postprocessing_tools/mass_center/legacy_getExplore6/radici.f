      subroutine radici(op, degree, zeror, zeroi)
      !op : vettore coefficienti polinomio 
      !degree: grado polinomio
      !zeror: vettore della parte reale delle radici
      !zeroi: vettore della parte immaginaria delle radici

      !trova le radici del polinomio di grado n ('degree') calcolando gli autovalori di
            !
            !  P(x)= a_1 x^(n)+ a_2 x^(n-1) + ........+ a_(n+1) x^(0)
            !
      ! M=[
      !  0   0                  -a_(n+1)/a_1
      !  1   0                  -a_(n  )/a_1
      !  0   1                  -a_(n-1)/a_1
      !  0   0   1              -a_(n-2)/a_1
      !                         
      !                         
      !  0   0   0          1   -a_(2  )/a_1
      !    ]
      integer degree,i,errore
      DOUBLE PRECISION op(degree+1),zeror(degree),zeroi(degree),M(degree,degree)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: fv1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iv1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: z

      ALLOCATE(z(degree,degree),iv1(degree),fv1(degree))
      !compongo matrice
            M(:,:)=0d0
            do i=1,degree
                  M(i,degree)=-op(degree+2-i)/op(1)
                  if (i>1) then
                        M(i,i-1)=1d0
                  endif
            enddo 
!            if ( rank .eq. 0) then
!            print*,'M ='
!            do i=1,degree
!                  print*,(M(i,j),j=1,degree)
!            enddo 
!            endif
      !cerco autovalori (routine eispack)
            call rg_l(degree,degree,M,zeror,zeroi,int(0),z,iv1,fv1,errore)
!            if ( rank .eq. 0) then
!                  print*,'zeror=',(zeror(j),j=1,degree)
!                  print*,'zeroi=',(zeroi(j),j=1,degree)
!                  stop
!            endif

      DEALLOCATE(z,iv1,fv1)
      end
