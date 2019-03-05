        PROGRAM evptn
!-------------------------------------------------------------------
! Use normalized density of sites to calculate the partion function 
! and then calculate the temperature-dependent vac formation 
! energies.
! <E_v> = c_A*<\epsilon_Av> + c_B*<\epsilon_Bv> + E^N/N
!                                                         Hu-Bin Luo
!-------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER :: i,j,k
        INTEGER :: iT,iTst,iTstp,iTend
        INTEGER :: na(2),ne(2)
        REAL(8),PARAMETER :: kb=0.0000063336231759d0   ! in Ry
        REAL(8) :: conc1,conc2
        REAL(8) :: ptn1,ptn2,eperf,eref1,eref2,eiv1,eiv2,eva,eva1,eva2
        CHARACTER(20) :: fnm(2)
        REAL(8) :: de(2),esub(2)
        REAL(8),ALLOCATABLE :: gga(:,:),ggb(:,:)
        REAL(8),ALLOCATABLE :: ffa(:),ffb(:)
!
!        WRITE(*,*) "The filenames of which where D(E)s are stored..."
        READ(*,*) (fnm(i),i=1,2)
        READ(*,*) (esub(i),i=1,2)
        READ(*,*) (na(i),i=1,2)
!        WRITE(*,*) "How many points of D(E) in the files..."
        READ(*,*) (ne(i),i=1,2)
!        WRITE(*,*) "Input the energy of perfect alloy (per site)..."
        READ(*,*) eperf
!        WRITE(*,*) "Tempertures: start, step and end..."
        READ(*,*) iTst,iTstp,iTend
!
        ALLOCATE(gga(ne(1),2),ggb(ne(2),2))
        ALLOCATE(ffa(ne(1)),ffb(ne(2)))
!
        DO i=1,2
           OPEN(8,file=TRIM(fnm(i)),status="old")
           DO j=1,ne(i)
              IF(i==1) READ(8,*) (gga(j,k),k=1,2)
              IF(i==2) READ(8,*) (ggb(j,k),k=1,2)
           END DO
           CLOSE(8)
        END DO
        conc1=1.d0*na(1)/(na(1)+na(2)); conc2=1.d0*na(2)/(na(1)+na(2))
!        gga(:,1)=gga(:,1)+esub(1)
!        ggb(:,1)=ggb(:,1)+esub(2)
!
! do integrals
!
        OPEN(9,file="evptn.prn",status="replace")
   !find out the eref to enhance the exp(-eiv/(kb*T))
        eref1=0.d0; eref2=0.d0
        DO i=1,ne(1)
           IF(gga(i,2)>1.d-3) THEN
              eref1=gga(i,1)
              EXIT
           END IF
        END DO
        DO j=1,ne(2)
           IF(ggb(j,2)>1.d-3) THEN
              eref2=ggb(j,1)
              EXIT
           END IF
        END DO
   !
        DO i=1,2
           IF(i==1) de(i)=(gga(ne(i),1)-gga(1,1))/(ne(i)-1)
           IF(i==2) de(i)=(ggb(ne(i),1)-ggb(1,1))/(ne(i)-1)
        END DO
   !calculate energies
        DO it=itst,itend,itstp
           DO i=1,ne(1)
                 eiv1=gga(i,1)+esub(1)
                 ffa(i)=gga(i,2)*eiv1*exp(-1.d0*(gga(i,1)-eref1)/(kb*it))
           END DO
           CALL gensim(ffa,de(1),ne(1),eva1)
           DO j=1,ne(2)
                 eiv2=ggb(j,1)+esub(2)
                 ffb(j)=ggb(j,2)*eiv2*exp(-1.d0*(ggb(j,1)-eref2)/(kb*it))
           END DO
           CALL gensim(ffb,de(2),ne(2),eva2)
   !calculate partion function
           DO i=1,ne(1)
                 ffa(i)=gga(i,2)*exp(-1.d0*(gga(i,1)-eref1)/(kb*it))
           END DO
           CALL gensim(ffa,de(1),ne(1),ptn1)
           DO j=1,ne(2)
                 ffb(j)=ggb(j,2)*exp(-1.d0*(ggb(j,1)-eref2)/(kb*it))
           END DO
           CALL gensim(ffb,de(2),ne(2),ptn2)
   !normalize the energies and print out
           eva1=eva1/ptn1; eva2=eva2/ptn2
           eva=conc1*eva1+conc2*eva2+eperf        
           WRITE(9,"(1X,I4,E16.6)") it,eva
        END DO
!
        CLOSE(9)
        DEALLOCATE(gga,ggb); DEALLOCATE(ffa,ffb)
!
        END PROGRAM       

!------------------------------------------------------------------------------

      SUBROUTINE gensim(f,d,n,fi)
!     ***************************************************************
!     *                                                             *
!     *    Integrates F from X(1) to X(n) where X is Louck's mesh.  *
!     *                                                             *
!     ***************************************************************
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(*) :: f
      REAL(KIND=8) :: d, fi, ff, corr
      INTEGER :: n, n1, j, n2
!
      IF(n.GE.3) THEN
         IF(MOD(n,2).EQ.0) THEN
           n1=n-1
           ff=-f(n-2)+8.d0*f(n-1)+5.d0*f(n)
           corr=ff*d/12.d0
         ELSE
           n1=n
           corr=0.d0
         ENDIF
         n2=n1-1
         ff=0.d0
         DO 20 j=2,n2,2
   20    ff=ff+(f(j-1)+4.d0*f(j)+f(j+1))
         fi=d*ff/3.d0+corr
      ELSEIF(n.EQ.2) THEN
         ff=5.d0*f(1)+8.d0*f(2)-f(3)
         fi=d*ff/12.d0
      ELSEIF(n.LE.1) THEN
         fi=0.d0
      ENDIF
!
      RETURN
      END
