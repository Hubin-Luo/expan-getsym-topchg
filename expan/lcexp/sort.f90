!--------------------------------------------------------------------
!
! subroutine to sort a vector elements in ascending order
!
        SUBROUTINE sort(vec,nelem)
        USE ctrl
        IMPLICIT NONE
        INTEGER :: nelem,i,j
        REAL(q) :: temp
        REAL(q),DIMENSION(nelem) :: vec
!
        DO i=1,nelem-1
           DO j=i+1,nelem
              IF(vec(i)>vec(j)) THEN
                temp=vec(i)
                vec(i)=vec(j)
                vec(j)=temp
              END IF
           END DO
        END DO
!
        END SUBROUTINE sort
