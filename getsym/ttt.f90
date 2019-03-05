        program ttt
        implicit none
        integer :: i, j, k
        real(8) :: a(3), b(3), c(3), ref(3), vmod(12)
        real(8) :: pi, quot
        
        READ(*,*)
        READ(*,*)
        DO i=1, 12
           READ(*,*)
        END DO
        READ(*,*)
        DO i=1, 12
           READ(*,*) (a(j), j=1,3), (b(j), j=1,3)
           vmod(i)=SQRT(DOT_PRODUCT(b,b))
        END DO
        quot=MAXVAL(vmod)
        vmod=vmod/quot*0.8
        WRITE(*,'(1X, 12F6.2)') (vmod(i), i=1, 12)

        end program
