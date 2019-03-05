        PROGRAM dmisym
!
        USE symmat
        IMPLICIT NONE
!
        INTEGER :: i, j, k, m, nbd, itarg
        INTEGER :: nsitein         ! number of input site
        REAL(8) :: modvec
        REAL(8) :: sympos(3),tarpos(3),vectmp(3)
        TYPE :: data
           INTEGER :: isite, nvec
           REAL(8) :: coor(3)
           REAL(8), ALLOCATABLE :: vecrel(:,:), dmi(:,:)
           REAL(8), ALLOCATABLE :: exch(:)
        END TYPE
        TYPE(data), ALLOCATABLE :: posin(:), posout(:)
!
        CALL getsym(1.d-5)
        ALLOCATE(posout(natom)) 
!
!read in number of central sites calculated from DFT
!
        READ(*,*) nsitein
        ALLOCATE(posin(nsitein))
!
!read in the site number and vector nubmer
! 
        DO i=1, nsitein
           READ(*,*) posin(i)%isite, posin(i)%nvec
           ALLOCATE(posin(i)%vecrel(3, posin(i)%nvec))
           ALLOCATE(posin(i)%dmi(3, posin(i)%nvec))
           ALLOCATE(posin(i)%exch(posin(i)%nvec))
           posin(i)%coor=pos(:, posin(i)%isite)
           DO j=1, posin(i)%nvec
              READ(*,*) (posin(i)%vecrel(k,j), k=1,3), (posin(i)%dmi(k,j), k=1,3),&
                       & posin(i)%exch(j)
              posin(i)%vecrel(:,j)=posin(i)%vecrel(:,j)-posin(i)%coor
           END DO
        END DO
!
!operate the central site and the vectors to map to all other equivalent
!positions
!
   !
   !operate the central site and check which position we get
   !
        DO i=1, natom
           posout(i)%coor=pos(:,i)
        END DO
        DO i=1, nsitein
           DO j=1, nop
              sympos=MATMUL(posin(i)%coor, mrot(:,:,j)) + vtran(:,j)
            !
            !translate if sympos is out of the cell
            !
              DO k=1, 3
                 IF (sympos(k)>1.d0) sympos(k)=sympos(k)-INT(sympos(k))
                 IF (sympos(k)<0.d0) sympos(k)=sympos(k)-INT(sympos(k))+1
              END DO
              Do k=1, natom
                 vectmp=sympos - posout(k)%coor
                 modvec=SQRT(DOT_PRODUCT(vectmp, vectmp))
                 IF (modvec > 1.d-5) CYCLE
              !
              !check if the position is already found by previous
              !operations. If not, put the vecrel and dmi in posout
              !
                 IF (.NOT. ALLOCATED(posout(k)%vecrel)) THEN 
                    ALLOCATE(posout(k)%vecrel(3, posin(i)%nvec))
                    ALLOCATE(posout(k)%dmi(3, posin(i)%nvec))
                    ALLOCATE(posout(k)%exch(posin(i)%nvec))
                    posout(k)%nvec=posin(i)%nvec
                    posout(k)%exch=posin(i)%exch
                    DO m=1, posout(k)%nvec
                       posout(k)%vecrel(:, m)=MATMUL(posin(i)%vecrel(:, m), mrot(:,:,j))
                       posout(k)%dmi(:, m)=MATMUL(posin(i)%dmi(:, m), mrot(:,:,j)) 
                    END DO
                    EXIT
                 END IF
              END DO
           END DO
        END DO
!
!write each position-related information
!
        OPEN(10, file='SITE-INFO', status='replace')
        DO i=1, natom
      !
      !check if posout(i)%vecrel is allocated. if so, print.
      !If not print zero information
      !
           IF (ALLOCATED(posout(i)%vecrel)) THEN
              WRITE(10, '(1X,I2,I4,3F15.9)') i, posout(i)%nvec, (posout(i)%coor(k),k=1,3)
              DO j=1, posout(i)%nvec
                !----------------------------------------------------
                ! find the site number of the target atom. If outside
                ! the cell, then translate into it
                !----------------------------------------------------
                 tarpos=posout(i)%coor+posout(i)%vecrel(:,j)
                 DO m=1, 3
                    IF(tarpos(m)>1.d0) tarpos(m)=tarpos(m)-INT(tarpos(m))
                    IF(tarpos(m)<0.d0) tarpos(m)=tarpos(m)-INT(tarpos(m))+1
                 END DO
                 DO m=1, natom
                    vectmp=tarpos-posout(m)%coor
                    modvec=SQRT(DOT_PRODUCT(vectmp, vectmp))
                    IF(modvec<1.d-5) THEN
                      itarg=m
                      EXIT
                    END IF
                 END DO
                 WRITE(10, '(3X, 2I4, 3F15.9, 4X, 3F15.9, 4X, F10.6)') &
               & i, itarg,                                             &
               & (posout(i)%vecrel(k, j), k=1,3),                      &
               & (posout(i)%dmi(k, j), k=1,3),                         &
               &  posout(i)%exch(j)
              END DO
           ELSE
              WRITE(10, '(1X, 2I4, 3F15.9)') i, 1, (posout(i)%coor(k),k=1,3)
              WRITE(10, '(1X, 3F15.9,4X,3F15.9,4X,F10.6)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
           END IF
        END DO
!
        END PROGRAM
