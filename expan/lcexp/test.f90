        program test
        implicit none
        integer:: i,ir
        real :: r 
        call init_seed()
        do i=1,10
          call random_number(r)
          ir=r*100+1
          write(*,*) ir
        enddo
        end program

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

