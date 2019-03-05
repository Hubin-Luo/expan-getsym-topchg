        PROGRAM evmc
!*********************************************************
! calculate the Evac using Monte Carlo simulation--H.B.L
!*********************************************************
        IMPLICIT NONE
        INTEGER,PARAMETER :: q=SELECTED_REAL_KIND(10)
        REAL(q),PARAMETER :: kb=0.0000063336231759_q   ! in Ry
        INTEGER :: i,j,k,m,it,ntot,isit
        INTEGER :: ncmp,ntst,itstp,ntend,nmc
        REAL(q) :: eperf,eiv,eiv1,eiv2,eva,ptn,etmp,blt,rand
        INTEGER,ALLOCATABLE :: natom(:)
        REAL(q),ALLOCATABLE :: esub(:)
        CHARACTER(50),ALLOCATABLE :: fnm(:)
        CHARACTER(1) :: opr
        TYPE :: spcsinfo 
          INTEGER,ALLOCATABLE :: site(:)
          REAL(q),ALLOCATABLE :: eng(:)
        END TYPE
        TYPE(spcsinfo),ALLOCATABLE :: spc(:)
!
        READ(*,*) ncmp
        ALLOCATE(fnm(ncmp))
        ALLOCATE(esub(ncmp)); ALLOCATE(natom(ncmp))
        READ(*,*) (fnm(i),i=1,ncmp)
        READ(*,*) (esub(i),i=1,ncmp)
        READ(*,*) (natom(i),i=1,ncmp)
        READ(*,*) eperf
        READ(*,*) opr,opr
        READ(*,*) ntst,itstp,ntend
        READ(*,*) nmc
        ntot=0
        DO i=1,ncmp
           ntot=ntot+natom(i)
        END DO
        ALLOCATE(spc(ncmp))
        DO i=1,ncmp
           ALLOCATE(spc(i)%site(natom(i)))
           ALLOCATE(spc(i)%eng(natom(i)))
        END DO
!
! read in the energies
!
        DO i=1,ncmp
           OPEN(7,file=TRIM(fnm(i)),status='old')
           DO j=1,6
              READ(7,*)
           END DO
           DO k=1,natom(i)
              READ(7,*) spc(i)%site(k),spc(i)%eng(k)
!              spc(i)%eng(k)=spc(i)%eng(k)+esub(i)
           END DO
           CLOSE(7)
        END DO
!
! opr=D: calculate directly. opr=M: Monte Carlo
!
        IF(opr=='D') THEN
           DO it=ntst,ntend,itstp
              eva=0._q; ptn=0._q
              DO i=1,natom(1)
                 DO j=1,natom(2)
                    eiv1=natom(1)*spc(1)%eng(i)+ &
                   &    natom(2)*spc(2)%eng(j)
                    eiv2=natom(1)*esub(1)+natom(2)*esub(2)
                    eiv1=eiv1-(eperf-INT(eperf))*(ntot-1)
                    eiv2=eiv2-INT(eperf)*(ntot-1)
                    eiv=eiv1+eiv2
                    eva=eva+eiv*exp(-1._q*eiv/(kb*it))
                    ptn=ptn+exp(-1._q*eva/(kb*it))
        write(*,*) "eiv: ", eiv
        write(*,*) "exp: ", eiv*exp(-1._q*eiv/(kb*it))
        write(*,*) "eva: ", eva
        write(*,*) "ptn: ", ptn
        stop
                 END DO
              END DO
              write(*,*) eva
              stop
              WRITE(*,"(1X,I4,2X,F11.6)") it,eva/ptn
           END DO
        ELSE IF(opr=='M') THEN
           CALL init_seed()
           DO it=ntst,ntend,itstp
              eva=0._q
              DO i=1,nmc
                 eiv1=0._q;eiv2=0._q;eiv=0._q
                 DO j=1,ncmp
                    CALL random_number(rand)
                    isit=INT(natom(j)*rand)+1
                    eiv1=eiv1+natom(j)*spc(j)%eng(isit)
                    eiv2=eiv2+natom(j)*esub(j)
                 END DO
                 eiv1=eiv1-(eperf-INT(eperf))*(ntot-1)
                 eiv2=eiv2-INT(eperf)*(ntot-1)
                 eiv =eiv1+eiv2
                 IF(i==1) THEN
                    etmp=eiv
                    eva=etmp/nmc
                    CYCLE
                 ELSE IF(eiv<etmp) THEN
                    eva=eva+eiv/nmc
                    etmp=eiv
                 ELSE IF(eiv>etmp) THEN
                    blt=exp((etmp-eiv)/(kb*it))
                    CALL random_number(rand)
                    IF(rand<blt) THEN
                       eva=eva+eiv/nmc
                       etmp=eiv
                    ELSE
                       eva=eva+etmp/nmc
                    END IF
                 END IF
              END DO
   !write the Monte Carlo results
   !           WRITE(*,"(1X,I4,2X,F10.6)") it,eva
               write(*,*) it,eva
           END DO
        ELSE
           WRITE(*,*) "No that operation"
           stop
        END IF
!
! deallocate the arrays
!
        DO i=1,ncmp
           DEALLOCATE(spc(i)%site)
           DEALLOCATE(spc(i)%eng)
        END DO
        DEALLOCATE(spc);  DEALLOCATE(fnm)
        DEALLOCATE(esub); DEALLOCATE(natom)
!
        END PROGRAM
!
!--------------------------------------------------
! init the seed of random number generator
!
        SUBROUTINE init_seed()
!
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL random_seed(size=n)
        ALLOCATE(seed(n))
!
        CALL system_clock(count=clock)
!
        seed=clock+37*(/(i-1,i=1,n)/)
        CALL random_seed(put = seed)

        DEALLOCATE(seed)
        END SUBROUTINE init_seed
