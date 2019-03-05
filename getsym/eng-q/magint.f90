        PROGRAM magint
!
! calculate magnetic interactions with a series of q state
!
        IMPLICIT NONE
        INTEGER :: i, j, k, iq, itmp
        INTEGER :: natom, nqp
        REAL(8) :: recip(3), qp(3), zz(3), xx(3)   !zz is the view axis
        REAL(8) :: spinorg(3),spintar(3), sptmp(3)
        REAL(8) :: lat(3,3)
        REAL(8) :: phase, modre, modzz, modxx
        REAL(8) :: qdr, edmitot, eexctot
        CHARACTER(1):: sprltyp, chrl
        CHARACTER(2):: fmt
        TYPE :: data
          INTEGER :: nvec
          REAL(8) :: coor(3)
          REAL(8), ALLOCATABLE :: edmi(:), eexc(:), exch(:)
          REAL(8), ALLOCATABLE :: vecrel(:,:), dmi(:,:)
        END TYPE
        TYPE(data), ALLOCATABLE :: pos(:)
!
        READ(*,*) natom
        READ(*,*) (recip(i), i=1,3) 
        READ(*,*) (zz(i), i=1,3)
        READ(*,*) nqp
        READ(*,*) sprltyp, phase
        READ(*,*) chrl
        DO i=1,3
           READ(*,*) (lat(i,j), j=1,3)
        END DO
        ALLOCATE(pos(natom))
        OPEN(8,file='SITE-INFO',status='old')
        DO i=1, natom
           READ(8,*) pos(i)%nvec, pos(i)%nvec, (pos(i)%coor(j),j=1,3)
           ALLOCATE(pos(i)%vecrel(3, pos(i)%nvec))
           ALLOCATE(pos(i)%dmi(3, pos(i)%nvec))
           ALLOCATE(pos(i)%exch(pos(i)%nvec))
           DO j=1, pos(i)%nvec
              READ(8,*) itmp, itmp, (pos(i)%vecrel(k,j), k=1,3), &
                       &(pos(i)%dmi(k,j), k=1,3),                &
                       & pos(i)%exch(j)
           END DO
           ALLOCATE(pos(i)%edmi(pos(i)%nvec))
           ALLOCATE(pos(i)%eexc(pos(i)%nvec))
        END DO
        CLOSE(8)
! 
! recip should be vertical to viewaxis. Find the axis recip crossp
! viewaxis
!
        IF (DOT_PRODUCT(recip,zz)<1.d-5) THEN
           CALL crossp(recip,zz,xx)
        ELSE
           WRITE(*,*) "ERROR: view axis is not vertical to q axis"
           STOP
        END IF
!
! calculate the magnetic interaction from site to site
!
        OPEN(9,file='ENG-INFO', status='replace')
        DO iq=0, nqp
          WRITE(9,'(1X, "iq=", I4)') iq
          edmitot=0.d0; eexctot=0.d0
          DO i=1, natom
          !
          !calculate the spin direction of coor
          !
             qdr=DOT_PRODUCT(recip/nqp*iq, pos(i)%coor)
             CALL spdir(qdr,phase,recip,xx,zz,sprltyp,chrl,spinorg)
             DO j=1, pos(i)%nvec
             !
             !calculate the spin directioan of vecrel
             !
                qdr=DOT_PRODUCT(recip/nqp*iq, pos(i)%coor+pos(i)%vecrel(:,j))
                CALL spdir(qdr,phase,recip,xx,zz,sprltyp,chrl,spintar)
             !
             !calculate the dmi energy first
             ! 
                CALL crossp(spinorg, spintar, sptmp)
                pos(i)%edmi(j)=-DOT_PRODUCT(pos(i)%dmi(:,j), sptmp)
             !
             !calculate the exchange energy then
             !
                pos(i)%eexc(j)=-pos(i)%exch(j)*13.6057*DOT_PRODUCT(spinorg, spintar)
             END DO
          !
          !write the interaction energies around the position
          !
             WRITE(fmt, '(I2)') pos(i)%nvec
             WRITE(9,'(1X, I4, 1X, "dmi and exch around", I4)') pos(i)%nvec, i
             WRITE(9,'(1X, '//TRIM(fmt)//'F10.4)') (pos(i)%edmi(k), k=1,pos(i)%nvec)
             WRITE(9,'(1X, '//TRIM(fmt)//'F10.4)') (pos(i)%eexc(k), k=1,pos(i)%nvec)
          !
          !calculate the total energies around the position
          !
             DO k=1,pos(i)%nvec
                edmitot=edmitot+pos(i)%edmi(k)
                eexctot=eexctot+pos(i)%eexc(k)
             END DO
          END DO
          WRITE(9,*)
          WRITE(9,'(8X, "Total magnetic interaction(meV):", F12.4)') (edmitot+eexctot)/2
          WRITE(9,*)
        END DO
        CLOSE(9)
!
        END PROGRAM
!
!==================================================================================
!
        SUBROUTINE spdir(qdr, phase, recip, xx, zz, sprltyp, chrl, spindir)
        IMPLICIT NONE
        INTEGER :: i, j, k
        REAL(8) :: recip(3), xx(3), zz(3), spindir(3)
        REAL(8) :: pi, phase, modre, modxx, modzz, qdr
        CHARACTER(1):: sprltyp, chrl
!
        pi=ACOS(-1.d0)
        modre=SQRT(DOT_PRODUCT(recip,recip))
        modzz=SQRT(DOT_PRODUCT(zz,zz))
        modxx=SQRT(DOT_PRODUCT(xx,xx))
!
        IF(sprltyp=="S") THEN
           IF(chrl=="R") THEN
              spindir=COS(2*pi*qdr+phase)*xx/modxx - &
                     &SIN(2*pi*qdr+phase)*zz/modzz
           ELSE IF(chrl=="L") THEN
              spindir=COS(2*pi*qdr+phase)*xx/modxx + &
                     &SIN(2*pi*qdr+phase)*zz/modzz
           ELSE
              WRITE(*,*) "No such kind of chirality"
              STOP
           END IF
        ELSE IF(sprltyp=="C") THEN
           IF(chrl=="R") THEN
              spindir=COS(2*pi*qdr+phase)*recip/modre - &
                     &SIN(2*pi*qdr+phase)*zz/modzz
           ELSE IF(chrl=="L") THEN
              spindir=COS(2*pi*qdr+phase)*recip/modre + &
                     &SIN(2*pi*qdr+phase)*zz/modzz
           ELSE
              WRITE(*,*) "No such kind of charity"
              STOP
           END IF
        ELSE
           WRITE(*,*) "No such kind of spiral type"
           STOP
        END IF 
!
        END SUBROUTINE              
