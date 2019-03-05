!==========================================================
!Cross product of two vectors                     H. B. Luo
!==========================================================
        SUBROUTINE crossp(vecfst, vecsec, vecout)
!
        IMPLICIT NONE
        INTEGER :: i, j, k
        REAL(8), INTENT(IN) :: vecfst(3), vecsec(3)
        REAL(8), INTENT(OUT):: vecout(3)
!
        vecout(1)=vecfst(2)*vecsec(3)-vecfst(3)*vecsec(2)
        vecout(2)=vecfst(3)*vecsec(1)-vecfst(1)*vecsec(3)
        vecout(3)=vecfst(1)*vecsec(2)-vecfst(2)*vecsec(1)
!
        END SUBROUTINE
