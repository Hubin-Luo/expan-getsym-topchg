        SUBROUTINE ineqvlt()
!
! tick out the inequivalent figures from sclst,pclst and tclst
! pos()%if() stores the positions of ineqvlt figs in [spt]clst
! Note that the numbers of columns of pos%fg's are set to be 2,
! 4 and 8, which are surely enough for binary because they 
! correspond to ternary alloy.
!
        USE ctrl; USE clsts; USE cellpos; 
        INTEGER :: i,j,k,m,n
        INTEGER :: itmp(50)
!
    !Single-point
        is=0; itmp=0
        DO i=1,ns
           DO j=1,i
              IF(ABS(sclst(i,4)-sclst(j,4))<thd) THEN
                 EXIT
              END IF
           END DO
           IF(i==j) THEN
              is=is+1; itmp(is)=i
           END IF
        END DO
        DO i=1,npos
           ALLOCATE(pos(i)%if1(is))
           ALLOCATE(pos(i)%fg1(is,2))
           pos(i)%fg1(:,:)=0._q
        END DO
        FORALL(i=1:npos)
           pos(i)%if1(:)=itmp(1:is)
        END FORALL
        WRITE(o1,*) "----inequivlt single-points----"
        DO i=1,is
           WRITE(o1,'(2X,I2,3F8.4,1X,F8.4)') i,(sclst(itmp(i),j),j=1,4)
        END DO
    !Pairs
        ip=0; itmp=0
        DO i=1,np
           DO j=1,i
              IF(MAXVAL(ABS(pclst(i,7:9)-pclst(j,7:9)))<thd) THEN
                 EXIT
              END IF
           END DO
           IF(i==j) THEN
              ip=ip+1; itmp(ip)=i
           END IF
        END DO
        DO i=1,npos
           ALLOCATE(pos(i)%if2(ip))
           ALLOCATE(pos(i)%fg2(ip,4))
           pos(i)%fg2(:,:)=0._q
        END DO
        FORALL(i=1:npos)
           pos(i)%if2(:)=itmp(1:ip)
        END FORALL
        WRITE(o1,*) "----inequivlt pairs----"
        DO i=1,ip
           WRITE(o1,130) i,(pclst(itmp(i),j),j=1,9)
           130 FORMAT(2X,I2,2(3F8.4,1X),2F8.4,1X,F8.4)
        END DO
    !Triplets
        it=0; itmp=0
        DO i=1,nt
           DO j=1,i
              IF(MAXVAL(ABS(tclst(i,10:15)-tclst(j,10:15)))<thd) THEN
                 EXIT
              END IF
           END DO
           IF(i==j) THEN
              it=it+1; itmp(it)=i
           END IF
        END DO
        DO i=1,npos
           ALLOCATE(pos(i)%if3(it))
           ALLOCATE(pos(i)%fg3(it,8))
           pos(i)%fg3(:,:)=0._q
        END DO
        FORALL(i=1:npos)
           pos(i)%if3(:)=itmp(1:it)
        END FORALL
        WRITE(o1,*) "----inequivlt triplets----"
        DO i=1,it
           WRITE(o1,140) i,(tclst(itmp(i),j),j=1,15)
           140 FORMAT(2X,I2,3(3F8.4,1X),3F8.4,1X,3F8.4)
        END DO
!
!get the no. of ECIs for fitting,plus 1 means constant term
!
        IF(ncomp==2) nterm=is+ip+it+1
        IF(ncomp==3) nterm=2*is+4*ip+8*it+1
        WRITE(o2,*)
        WRITE(o2,"(1X,'Number of terms(plus intercept): ',I4)") nterm
        WRITE(o2,*)
!
        END SUBROUTINE
