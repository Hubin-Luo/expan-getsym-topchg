        SUBROUTINE sor(a,b,n,w,k,x)
!------------------------------------------
! subroutine to iterate using SOR algorithm
! which requires a(i,i) be a non-zero value
! to solve equations Ax=b
!------------------------------------------
        USE ctrl; 
        IMPLICIT NONE
        INTEGER :: i,j,m,n
        INTEGER,INTENT(out):: k
        REAL(q),INTENT(in) :: w
        REAL(q),INTENT(in) :: a(n,n)
        REAL(q),INTENT(in) :: b(n)
        REAL(q),INTENT(out):: x(n)
        REAL(q) :: dx(n)
        REAL(q) :: sum
!
        k=1; x=0._q; dx=0._q
        DO WHILE (k<20000)
           DO i=1,n
              sum=0._q
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
