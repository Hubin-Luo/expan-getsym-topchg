        MODULE cellpos
!------------------------------------------------------------------------
! if1-3 stores the positions of ineqvlt figs in [spt]clst (see
! ineqvlt.f90). fg1-3 stores the correlation functions (see correl.f90)
!------------------------------------------------------------------------
        USE ctrl
        INTEGER :: npos,is,ip,it
        TYPE :: position
           INTEGER :: iatom
           CHARACTER(LEN=1) :: atype
           INTEGER :: iq
           REAL(q),DIMENSION(3) :: coor
           INTEGER,ALLOCATABLE :: if1(:),if2(:),if3(:)
           REAL(q),ALLOCATABLE :: fg1(:,:),fg2(:,:),fg3(:,:) 
        END TYPE
        TYPE(position),ALLOCATABLE :: pos(:)
!
        CONTAINS
!
! read in the positions in the supercell
!
        SUBROUTINE getpos()
!
        IMPLICIT NONE
        INTEGER :: i,j
        CHARACTER(LEN=150) :: line
!
        DO WHILE(.TRUE.)
           READ(i1,'(A)') line
           IF (line(12:17)=='isc at') EXIT
        END DO
  !count the number of position
        npos=0
        DO WHILE(.TRUE.)
           READ(i1,'(A)') line
           IF (line(24:41)=='                  ') THEN
              EXIT
           ELSE
              npos=npos+1
           END IF
        END DO
        ALLOCATE(pos(npos))
  !backspace to read the positions
        DO 10 i=1,npos+1
     10    BACKSPACE(UNIT=i1)
        DO 20 i=1,npos
     20    READ(i1,*) pos(i)%iatom, pos(i)%atype, pos(i)%iq, &
      &              (pos(i)%coor(j),j=1,3)
  !feed back the positions
  !     WRITE(o1,*) "           isc at  IQ    QXSC     QYSC     QZSC"
  !     DO i=1,npos
  !        WRITE(o1,100) pos(i)%iatom, pos(i)%atype, pos(i)%iq,&
  !   &                 (pos(i)%coor(j),j=1,3)
  !        100  FORMAT(10X,I4,1X,A1,1X,I4,3(1X,F9.6))
  !     END DO
!
! It is convenient to assign the occupation number here to iatom variable
!
        DO i=1,npos
           IF(pos(i)%atype=="A") pos(i)%iatom=1
           IF(pos(i)%atype=="B") pos(i)%iatom=-1
           IF(pos(i)%atype=="C") pos(i)%iatom=0
        END DO
!        DO i=1,npos
!           WRITE(6,*) pos(i)%iatom
!        END DO
!
        END SUBROUTINE
!
        END MODULE cellpos
