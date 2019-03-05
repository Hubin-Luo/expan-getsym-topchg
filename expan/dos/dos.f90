        PROGRAM dos
!-----------------------------------------------------------
! calculate the dos using Gaussion function 
! It reads the file with surrogate energies output from 
! expan.exe calcualtion                           Hu-Bin Luo
!-----------------------------------------------------------
        IMPLICIT NONE
        INTEGER :: i,ie,nein,neout,itmp
        REAL(8) :: pi
        REAL(8),ALLOCATABLE :: eng(:)
        REAL(8) :: sigma,emin,emax,de,eout,ediff
        REAL(8) :: tdos
        CHARACTER(20) :: fnm
!
        pi=ACOS(-1.d0)
        WRITE(*,*) "input the filename where energies stored..."
        READ(*,*) fnm
        WRITE(*,*) "numbers of input and output energies..."
        READ(*,*) nein, neout
        WRITE(*,*) "input the smearing in Gaussion function..."
        READ(*,*) sigma
        WRITE(*,*) "input emin, emax..."
        READ(*,*) emin,emax
!
        ALLOCATE(eng(nein))
!
        OPEN(8,file=TRIM(fnm),status="old")
        OPEN(9,file="dos.txt",status="replace")
        DO i=1,6
           READ(8,*)
        END DO
        DO i=1,nein
           READ(8,*) itmp,eng(i)
        END DO
        CLOSE(8)
!
        de=(emax-emin)/(neout-1)
        DO ie=1,neout
           eout=emin+(ie-1)*de
           tdos=0.d0
           DO i=1,nein
              ediff=eng(i)-eout
              tdos=tdos+exp(-ediff*ediff/(2*sigma*sigma)) &
             &        /SQRT(2*pi*sigma*sigma)
           END DO
           tdos=tdos/nein
           WRITE(9,*) eout, tdos
        END DO
        CLOSE(9)
        DEALLOCATE(eng)
!
        END PROGRAM
