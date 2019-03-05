        program test
        use symmat
        implicit none
        integer :: i,j,k
        real(8) :: ptmp(3), vtmp(3), symv(3), sympos(3)
        real(8) :: tmp(3), vectmp1(3), vectmp2(3)
        real(8) :: bndmod

        ptmp(1)=0.05670 
        ptmp(2)=0.05670 
        ptmp(3)=0.05670 
        sympos=0.d0
        CALL getsym(1.0d-5)
        do i=1, nop
           sympos=MATMUL(ptmp, mrot(:,:,i)) + vtran(:,i)
           do k=1,3
              if (sympos(k)>1.d0) sympos(k)=sympos(k)-INT(sympos(k))
              if (sympos(k)<0.d0) sympos(k)=sympos(k)-INT(sympos(k))+1
           end do  
              do j=1, natom
                   tmp=sympos-pos(:,j)
                   if (dot_product(tmp,tmp) < 1.d-4) then
                      write(*,*) "operation", i, "transformed to position", j
                   end if
                end do
        end do
        
        end program
           
