        SUBROUTINE dos(eng, sigma, tdos, nein, neout)
!-----------------------------------------------------------
! calculate the dos using Gaussion function 
!                                Hu-Bin Luo
!-----------------------------------------------------------
        IMPLICIT NONE
        INTEGER :: i,ie,nein,neout
        REAL(8) :: pi
        REAL(8) :: eng(nein)
        REAL(8) :: sigma,emin,emax,de,ediff
        REAL(8),ALLOCATABLE,INTENT(OUT) :: tdos(:,:)
!
        pi=ACOS(-1.d0)
        emin=MINVAL(eng)-(MAXVAL(eng)-MINVAL(eng))*0.2
        emax=MAXVAL(eng)+(MAXVAL(eng)-MINVAL(eng))*0.2
!
        ALLOCATE(tdos(neout,2))
!
        tdos(:,:)=0.d0
        de=(emax-emin)/(neout-1)
        DO ie=1,neout
           tdos(ie,1)=emin+(ie-1)*de
           DO i=1,nein
              ediff=eng(i)-tdos(ie,1)
              tdos(ie,2)=tdos(ie,2)+exp(-ediff*ediff/(2*sigma*sigma)) &
             &        /SQRT(2*pi*sigma*sigma)
           END DO
           tdos(ie,2)=tdos(ie,2)/nein
        END DO
        DEALLOCATE(tdos)
!
        END SUBROUTINE
