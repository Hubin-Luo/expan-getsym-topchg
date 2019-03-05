        SUBROUTINE sor(a,b,n,w,k,x)
!------------------------------------------
! subroutine to iterate using SOR algorithm
! which requires a(i,i) be a non-zero value
! to solve equations Ax=b
!------------------------------------------
        IMPLICIT NONE
        INTEGER :: i,j,m,n
        INTEGER,INTENT(out):: k
        REAL(8),INTENT(in) :: w
        REAL(8),INTENT(in) :: a(n,n)
        REAL(8),INTENT(in) :: b(n)
        REAL(8),INTENT(out):: x(n)
        REAL(8) :: dx(n)
        REAL(8) :: sum, hhd
!
        hhd=1.d-6
        k=1; x=0.d0; dx=0.d0
        DO WHILE (k<20000)
           DO i=1,n
              sum=0.d0
              DO j=1,n
                 sum=sum+a(i,j)*x(j)
              END DO
              dx(i)=w/a(i,i)*(b(i)-sum)
              x(i)=x(i)+dx(i)
           END DO
           IF(MAXVAL(ABS(dx(1:n)))<hhd) EXIT
           k=k+1
        END DO
!
        END SUBROUTINE          
