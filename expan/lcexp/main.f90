        PROGRAM EXPAN
        
!*********************************************************************
!Description:                                                        *
!This code is to calculate the correlation function around each site *
!and then fit the energies to calculate ECI. Then surrogate energies *
!can thus be calculated through expansion.                           *
!                                          ---------Hu-Bin Luo       *
!*********************************************************************
!
        USE ctrl; USE clsts; USE cellpos; 
        USE abin
!
        IMPLICIT NONE
        INTEGER :: i,j,k,m,n,iexch
        CHARACTER(1) :: opr,spc,constc
        CHARACTER(150):: line,chtmp
        REAL(q),ALLOCATABLE :: eci(:)
!
! read in the control inputs for calculation of correlation function
!
        READ(*,*) fnm, ncomp
        READ(*,*) chtmp, constc
        IF(constc=='Y') THEN
          READ(*,*) chtmp, iexch
          WRITE(*,*)
          WRITE(*,*) "Composition constant! vac site:", iexch
          WRITE(*,*)
        ELSE
          READ(*,*)
        END IF
        READ(*,*) (bs1(i),i=1,3)
        READ(*,*) (bs2(i),i=1,3)
        READ(*,*) (bs3(i),i=1,3)
        READ(*,*) (abc(i),i=1,3)
        READ(*,*) (nbs(j),j=1,3)
        READ(*,*) csize,slen
        READ(*,*) (plen(i),i=1,2)
        READ(*,*) (tlen(i),i=1,2)
!
! read in the control parameter what to do after calculation of cor.
! func.
!
        READ(*,*)
        READ(*,*) opr,opr
        READ(*,*)
!
! read in the ab inito energies or ECIs
!
!        IF(opr=='F') THEN
!           READ(*,*) ne, nesub
!           ALLOCATE(ist(ne),eng(ne))
!           DO i=1,ne
!              READ(*,*) ist(i),eng(i)
!              eng(i)=eng(i)-nesub
!           END DO
!        ELSE IF(opr=='E') THEN
!           DO WHILE(.TRUE.)
!              READ(*,'(A)') line
!              IF(line(1:3)=='---') EXIT
!           END DO
!           READ(*,*) spc,spc
!           ALLOCATE(eci(neci))
!           READ(*,*) (eci(i),i=1,neci)
!        END IF
              
!
        OPEN(i1,file=TRIM(fnm),status='old')
        OPEN(o1,file='POSN',status='replace')
        OPEN(o2,file='CRLF',status='replace')
        OPEN(o3,file='FTEX',status='replace')
!
! read in the positions and generate the clusters
!
        CALL getpos()
        CALL gclst()
!
! tick out the inequivalent figures
!
        CALL ineqvlt()
!
! run out the positions to get all the correlation functions
!
        WRITE(o2,*) "*********************************"
        WRITE(o2,*) "Calculated correlation functions:"
        WRITE(O2,*) "*********************************"
        IF(constc=='Y') THEN
           DO i=1,npos
              pos(iexch)%iatom=pos(i)%iatom
              CALL correl(i)
           END DO
        ELSE
           DO i=1,npos
              CALL correl(i) 
           END DO
        END IF
!
! Two cases for operation. One is fitting and another
! is to calculate surrogate energies
!
        IF(opr=='F') THEN 
  !
           READ(*,*) ne, nesub
           ALLOCATE(ist(ne),eng(ne))
           DO i=1,ne
              READ(*,*) ist(i),eng(i)
              eng(i)=eng(i)-nesub
           END DO
  !          
           CALL fit()
  !
        ELSE IF(opr=='E') THEN
  !
           DO WHILE(.TRUE.)
              READ(*,'(A)') line
              IF(line(1:3)=='---') EXIT
           END DO
           READ(*,*) spc,spc
  !
           ALLOCATE(eci(nterm))
           READ(*,*) (eci(i),i=1,nterm)
  !
           CALL engcal(spc,eci)
  !
        ELSE
           WRITE(*,*) "  Only correlation functions calculated  "
        END IF        
!
! close the units and deallocate the arrays
!
        DO i=1,npos
           DEALLOCATE(pos(i)%if1,pos(i)%if2,pos(i)%if3)
           DEALLOCATE(pos(i)%fg1,pos(i)%fg2,pos(i)%fg3)
        END DO
        DEALLOCATE(pos)
        CLOSE(i1)
        CLOSE(o1)
        CLOSE(o2)
        CLOSE(o3)
!
        END PROGRAM
