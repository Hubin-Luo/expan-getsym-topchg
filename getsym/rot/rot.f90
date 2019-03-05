        PROGRAM rot
        IMPLICIT NONE
        INTEGER :: i, j, k, m, nvec, nsite
        REAL(8) :: pi,  thetadeg
        REAL(8) :: vecin(3), vecout(3), dmiin(3), dmiout(3), axis(3), trv(3)
!
        OPEN(8,file='INPUTCAL', status='old')
        READ(8,*) nsite
        DO i=1, nsite
           READ(8,*) nvec, (axis(j), j=1,3), (trv(j), j=1,3), thetadeg
           WRITE(*,'(1X, I4, 2X, F7.2)') nvec, thetadeg
           DO j=1, nvec
              READ(8,*) (vecin(k), k=1,3), (dmiin(m), m=1,3)
              CALL rotardaxis(axis, trv, thetadeg, vecin, vecout, 0)
              CALL rotardaxis(axis, trv, thetadeg, dmiin, dmiout, 1)
              WRITE(*,'(1X, 3F12.6, 4X, 3F12.6)') (vecout(k), k=1,3), &
             &     (dmiout(m), m=1,3)
           END DO
        END DO
!
        END PROGRAM
