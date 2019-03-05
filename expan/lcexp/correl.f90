        SUBROUTINE correl(ipos)
!--------------------------------------------------------
! get the correlation functions around the specified site
! taking the multicomponent system into account. 
! Note that when print out the correlation functions, the 
! corresponding figures are distermined in ineqvlt.f90. See
! the outputs in "POSN"
!--------------------------------------------------------
        USE ctrl; USE cellpos; USE clsts;
        IMPLICIT NONE
        INTEGER :: ipos,i,j,k,ii
        REAL(q) :: stmp(ns,4),ptmp(np,9),ttmp(nt,15)
        REAL(q) :: op1,op2,op3
        REAL(q) :: ot1(2),ot2(2),ot3(2)
        CHARACTER(3) :: fm1,fm2,fm3
        REAL(q),EXTERNAL :: thet1,thet2
!
!transfer the clst to tmps
!
        stmp=sclst(1:ns,:)
        ptmp=pclst(1:np,:)
        ttmp=tclst(1:nt,:) 
!
!translate the position
!
        DO i=1,ns
           stmp(i,1:3)=stmp(i,1:3)+pos(ipos)%coor(:)
        END DO
        DO i=1,np
           ptmp(i,1:3)=ptmp(i,1:3)+pos(ipos)%coor(:)
           ptmp(i,4:6)=ptmp(i,4:6)+pos(ipos)%coor(:)
        END DO
        DO i=1,nt
           ttmp(i,1:3)=ttmp(i,1:3)+pos(ipos)%coor(:)
           ttmp(i,4:6)=ttmp(i,4:6)+pos(ipos)%coor(:)
           ttmp(i,7:9)=ttmp(i,7:9)+pos(ipos)%coor(:)
        END DO
!
!make sure all the surrounding atoms are in the cell
!
        FORALL(i=1:ns,j=1:3,stmp(i,j)<0._q)           &
         & stmp(i,j)=stmp(i,j)+ &
         & INT(1+ABS(stmp(i,j)/abc(j)))*abc(j)
        FORALL(i=1:ns,j=1:3,stmp(i,j)>=abc(j))        &
         & stmp(i,j)=stmp(i,j)- &
         & INT(ABS(stmp(i,j)/abc(j)))*abc(j)
        FORALL(i=1:np,j=1:6,ptmp(i,j)<0._q)           &
         & ptmp(i,j)=ptmp(i,j)+ &
         & INT(1+ABS(ptmp(i,j)/abc(j)))*abc(j-3*INT((j-0.5_q)/3))
        FORALL(i=1:np,j=1:6,ptmp(i,j)>=abc(j-3*INT((j-0.5_q)/3)))  &
         & ptmp(i,j)=ptmp(i,j)- &
         & INT(ABS(ptmp(i,j)/abc(j)))*abc(j-3*INT((j-0.5_q)/3))
        FORALL(i=1:nt,j=1:9,ttmp(i,j)<0._q)           &
         & ttmp(i,j)=ttmp(i,j)+ &
         & INT(1+ABS(ttmp(i,j)/abc(j)))*abc(j-3*INT((j-0.5_q)/3))
        FORALL(i=1:nt,j=1:9,ttmp(i,j)>=abc(j-3*INT((j-0.5_q)/3)))  &
         & ttmp(i,j)=ttmp(i,j)- &
         & INT(ABS(ttmp(i,j)/abc(j)))**abc(j-3*INT((j-0.5_q)/3))
    !DO i=1,np
    !   WRITE(6,110) (ptmp(i,j),j=1,9)
    !     110 FORMAT(2X,2(3F8.4,1X),2F8.4,1X,F8.4)
    !END DO
!
!----------------------------------------------------------------------
!find the occ. no. and calculate the correlation func.
!Note that the occ. nos. are stored in pos(:)%iatom variables, so it is
!necessary to compare the positions in [spt]tmp with pos(:)%coor. It is 
!also important to add up the corr. func. of equivalent figures.
!After the functions are sucessfully calculated and stored in pos(:)%fg
!, write it to printout file attached to unit o2
!----------------------------------------------------------------------
!
!write the number of inequivalent figs. to fm variable and thus can be
!used to print out the correlation functions
!
        fm1="??"; fm2="??"; fm3="??"
        WRITE(fm1,'(I2)') is
        WRITE(fm2,'(I2)') ip
        WRITE(fm3,'(I3)') it
!
!calculate correlation functions
!
        SELECT CASE(ncomp)
!
!binary--the cluster functions are simply constructed as A->1,B->-1
!
        CASE(2)
   !single-points
        DO i=1,is
           ii=pos(ipos)%if1(i)
           DO j=ii,ns
              IF(ABS(stmp(ii,4)-stmp(j,4))<thd) THEN
                 DO k=1,npos
                    IF(MAXVAL(ABS(stmp(j,1:3)-pos(k)%coor(:)))<thd) &
                    &  op1=pos(k)%iatom
                 END DO
                 pos(ipos)%fg1(i,1)=pos(ipos)%fg1(i,1)+op1
              END IF
           END DO
        END DO
   !pairs
        DO i=1,ip
           ii=pos(ipos)%if2(i)
           DO j=ii,np
              IF(MAXVAL(ABS(ptmp(ii,7:9)-ptmp(j,7:9)))<thd) THEN
                 DO k=1,npos
                    IF(MAXVAL(ABS(ptmp(j,1:3)-pos(k)%coor(:)))<thd) &
                    &   op1=pos(k)%iatom
                    IF(MAXVAL(ABS(ptmp(j,4:6)-pos(k)%coor(:)))<thd) &
                    &   op2=pos(k)%iatom
                 END DO
                 pos(ipos)%fg2(i,1)=pos(ipos)%fg2(i,1)+op1*op2
              END IF   
           END DO
        END DO
   !triplets
        DO i=1,it
           ii=pos(ipos)%if3(i)
           DO j=ii,nt
              IF(MAXVAL(ABS(ttmp(ii,10:15)-ttmp(j,10:15)))<thd) THEN
                 DO k=1,npos
                    IF(MAXVAL(ABS(ttmp(j,1:3)-pos(k)%coor(:)))<thd) &
                    &  op1=pos(k)%iatom
                    IF(MAXVAL(ABS(ttmp(j,4:6)-pos(k)%coor(:)))<thd) &
                    &  op2=pos(k)%iatom
                    IF(MAXVAL(ABS(ttmp(j,7:9)-pos(k)%coor(:)))<thd) &
                    &  op3=pos(k)%iatom
                 END DO
                 pos(ipos)%fg3(i,1)=pos(ipos)%fg3(i,1)+op1*op2*op3
              END IF
           END DO
        END DO
   !write the functions to printout file
        WRITE(o2,'(1X,I4,'//fm1//'(F7.2),1X,"+",    &
        & '//fm2//'(F7.2),1X,"+",'//fm3//'(F7.2))') &
        & ipos,                                     &
        & (pos(ipos)%fg1(i,1),i=1,is),              &
        & (pos(ipos)%fg2(i,1),i=1,ip),              &
        & (pos(ipos)%fg3(i,1),i=1,it)
!
!ternary--the cluster functions are constructed by Chebyshev polynomials
! 
        CASE(3)
   !single-points
        DO i=1,is
           ii=pos(ipos)%if1(i)
           DO j=ii,ns
              IF(ABS(stmp(ii,4)-stmp(j,4))<thd) THEN
                 DO k=1,npos
                    IF(MAXVAL(ABS(stmp(j,1:3)-pos(k)%coor(:)))<thd) THEN
                       ot1(1)=thet1(pos(k)%iatom)
                       ot1(2)=thet2(pos(k)%iatom)
                    END IF
                 END DO
                 pos(ipos)%fg1(i,1)=pos(ipos)%fg1(i,1)+ot1(1)
                 pos(ipos)%fg1(i,2)=pos(ipos)%fg1(i,2)+ot1(2)
              END IF
           END DO
        END DO
   !pairs
        DO i=1,ip
           ii=pos(ipos)%if2(i)
           DO j=ii,np
              IF(MAXVAL(ABS(ptmp(ii,7:9)-ptmp(j,7:9)))<thd) THEN
                 DO k=1,npos
                    IF(MAXVAL(ABS(ptmp(j,4:6)-pos(k)%coor(:)))<thd) THEN
                       ot1(1)=thet1(pos(k)%iatom)
                       ot1(2)=thet2(pos(k)%iatom)
                    END IF
                    IF(MAXVAL(ABS(ptmp(j,7:9)-pos(k)%coor(:)))<thd) THEN
                       ot2(1)=thet1(pos(k)%iatom)
                       ot2(2)=thet2(pos(k)%iatom)
                    END IF
                 END DO
                 pos(ipos)%fg2(i,1)=pos(ipos)%fg2(i,1)+ot1(1)*ot2(1)
                 pos(ipos)%fg2(i,2)=pos(ipos)%fg2(i,2)+ot1(1)*ot2(2)
                 pos(ipos)%fg2(i,3)=pos(ipos)%fg2(i,3)+ot1(2)*ot2(1)
                 pos(ipos)%fg2(i,4)=pos(ipos)%fg2(i,4)+ot1(2)*ot2(2)
              END IF
           END DO
        END DO
   !triplets
        DO i=1,it
           ii=pos(ipos)%if3(i)
           DO j=ii,nt
              IF(MAXVAL(ABS(ttmp(ii,10:15)-ttmp(j,10:15)))<thd) THEN
                 DO k=1,npos
                    IF(MAXVAL(ABS(ttmp(j,1:3)-pos(k)%coor(:)))<thd) THEN
                       ot1(1)=thet1(pos(k)%iatom)
                       ot1(2)=thet2(pos(k)%iatom)
                    END IF
                    IF(MAXVAL(ABS(ttmp(j,4:6)-pos(k)%coor(:)))<thd) THEN
                       ot2(1)=thet1(pos(k)%iatom)
                       ot2(2)=thet2(pos(k)%iatom)
                    END IF
                    IF(MAXVAL(ABS(ttmp(j,7:9)-pos(k)%coor(:)))<thd) THEN
                       ot3(1)=thet1(pos(k)%iatom)
                       ot3(2)=thet2(pos(k)%iatom)
                    END IF
                 END DO
                 pos(ipos)%fg3(i,1)=pos(ipos)%fg3(i,1)+ot1(1)*ot2(1)*ot3(1)
                 pos(ipos)%fg3(i,2)=pos(ipos)%fg3(i,2)+ot1(1)*ot2(1)*ot3(2)
                 pos(ipos)%fg3(i,3)=pos(ipos)%fg3(i,3)+ot1(1)*ot2(2)*ot3(1)
                 pos(ipos)%fg3(i,4)=pos(ipos)%fg3(i,4)+ot1(1)*ot2(2)*ot3(2)
                 pos(ipos)%fg3(i,5)=pos(ipos)%fg3(i,5)+ot1(2)*ot2(1)*ot3(1)
                 pos(ipos)%fg3(i,6)=pos(ipos)%fg3(i,6)+ot1(2)*ot2(1)*ot3(2)
                 pos(ipos)%fg3(i,7)=pos(ipos)%fg3(i,7)+ot1(2)*ot2(2)*ot3(1)
                 pos(ipos)%fg3(i,8)=pos(ipos)%fg3(i,8)+ot1(2)*ot2(2)*ot3(2)
              END IF
           END DO
        END DO
   !write the functions to printout file
        WRITE(o2,'(1X,I4,'//fm1//'(2F7.2),1X,"+",     &
        & '//fm2//'(4F7.2),1X,"+",'//fm3//'(8F7.2))') &
        & ipos,                                       &
        & ((pos(ipos)%fg1(i,j),j=1,2),i=1,is),        &
        & ((pos(ipos)%fg2(i,j),j=1,4),i=1,ip),        &
        & ((pos(ipos)%fg3(i,j),j=1,8),i=1,it)

!
        CASE(4:)
          WRITE(6,*) "not suported beyond ternary now"
          STOP
!
        END SELECT
!
        END SUBROUTINE
!
!the polynomials of ternary system
!
        FUNCTION thet1(ix)
        USE ctrl
        IMPLICIT NONE
        INTEGER :: ix
        REAL(q) :: thet1
        thet1=SQRT(1.5_q)*ix
        END FUNCTION
!
        FUNCTION thet2(ix)
        USE ctrl
        IMPLICIT NONE
        INTEGER :: ix
        REAL(q) :: thet2
        thet2=SQRT(2._q)*(1._q - 1.5_q*ix*ix)
        END FUNCTION
