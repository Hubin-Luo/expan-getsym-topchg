        SUBROUTINE engcal(spc,eci)
        USE ctrl; USE cellpos;
        USE clsts
        IMPLICIT NONE
        INTEGER :: i,j,k,m
        REAL(q) :: eci(nterm),fun(nterm)
        REAL(q) :: fen
        CHARACTER(1) :: spc
        CHARACTER(3) :: fmt
!
        WRITE(o3,*) "**********************"
        WRITE(o3,*) "The surrogate energies"
        WRITE(o3,*) "**********************"
        WRITE(o3,*) "Atom type: ", spc
        fmt="???"
        WRITE(fmt,'(I3)') nterm
        WRITE(o3,"(1X,'ECIs: ',"//fmt//"(E16.8)/)") &
              &(eci(i),i=1,nterm)
!
        fun(1)=1._q
!
        IF(ncomp==2) THEN
           DO i=1,npos
              IF(pos(i)%atype==spc) THEN
                 DO j=2,nterm
                    IF(j<=is+1) THEN
                       fun(j)=pos(i)%fg1(j-1,1)
                    ELSE IF(j<=is+ip+1) THEN
                       fun(j)=pos(i)%fg2(j-1-is,1)
                    ELSE
                       fun(j)=pos(i)%fg3(j-1-is-ip,1)
                    END IF
                 END DO
                 fen=DOT_PRODUCT(eci(:),fun(:))
                 WRITE(o3,"(1X,I4,1X,F12.6)") i,fen
              END IF
           END DO
        END IF
!
        IF(ncomp==3) THEN
           DO i=1,npos
              IF(pos(i)%atype==spc) THEN
                 DO j=2,nterm
                    IF(j<=2*is+1) THEN
                       k=INT( (j-1+1.5_q)/2 )
                       m=j-2*INT((j-1-0.5_q)/2)
                       fun(j)=pos(i)%fg1(k,m)
                    ELSE IF(j<=2*is+4*ip+1) THEN
                       k=INT( (j-1-2*is+1.5_q)/4 )
                       m=j-1-2*is-4*INT( (j-1-2*is-0.5_q)/4 )
                       fun(j)=pos(i)%fg2(k,m)
                    ELSE
                       k=INT( (j-1-2*is-4*ip+1.5_q)/8 )
                       m=j-1-2*is-4*ip-8*INT( (j-1-2*is-4*ip)/8 )
                       fun(j)=pos(i)%fg3(k,m)
                    END IF
                 END DO
                 fen=DOT_PRODUCT(eci(:),fun(:))
                 WRITE(o3,"(1X,I4,1X,F12.6)") i,fen
              END IF
           END DO
        END IF
!
        END SUBROUTINE
