        MODULE clsts
!
        USE ctrl
        INTEGER :: nv,ns,np,nt
        REAL(q) :: csize,slen
        INTEGER,DIMENSION(3) :: nbs=(/1,1,1/)
        REAL(q),DIMENSION(3) :: abc=(/1.,1.,1./)
        REAL(q),DIMENSION(2) :: plen, tlen
        REAL(q),DIMENSION(3) :: bs1,bs2,bs3
        REAL(q),DIMENSION(200,3) :: vcs
        REAL(q),DIMENSION(200,4) :: sclst
        REAL(q),DIMENSION(500,9) :: pclst
        REAL(q),DIMENSION(500,15):: tclst
!
        CONTAINS
!
! choose one position (0,0,0) as the center and find the inequilibrium
! clusters that satisfy the condition within the certain shell
!
        SUBROUTINE gclst()
!        
        USE ctrl; 
        IMPLICIT NONE
        INTEGER :: i,j,k,m
        REAL(q) :: dist1,dist2,dist3
        REAL(q),DIMENSION(3) :: vtp1,vtp2,vtp3
!        
        vcs(:,:)=0._q; pclst(:,:)=0._q; tclst(:,:)=0._q
!
  !generate vectors satisfied
        nv=0
        DO i=ABS(nbs(1)),-1*ABS(nbs(1)),-1
           DO j=ABS(nbs(2)),-1*ABS(nbs(2)),-1
              DO k=ABS(nbs(3)),-1*ABS(nbs(3)),-1
                 DO m=1,3
                    vtp1(m)=i*bs1(m)+j*bs2(m)+k*bs3(m)
                 END DO
                 dist1=SQRT(DOT_PRODUCT(vtp1,vtp1))
                 IF(dist1>thd .AND. dist1<=csize) THEN
                    nv=nv+1
                    vcs(nv,:)=vtp1(:)
                 END IF
              END DO
           END DO
        END DO
        WRITE(o1,*) "Local positions:"
        DO i=1,nv
           WRITE(o1,'(2X,3F8.4)') (vcs(i,j),j=1,3)
        END DO
  !single-point clst
        ns=0
        DO i=1,nv
           dist1=SQRT(DOT_PRODUCT(vcs(i,:),vcs(i,:)))
           IF(dist1<=slen) THEN
             ns=ns+1
             sclst(ns,1:3)=vcs(i,:)
             sclst(ns,4)=dist1
           END IF
        END DO
        WRITE(o1,*) "----Single-point----"
        DO i=1,ns
           WRITE(o1,'(2X,3F8.4,1X,F8.4)') (sclst(i,j),j=1,4)
        END DO
  !pairs
        np=0
        DO i=1,nv-1
           DO j=i+1,nv
              vtp1(:)=vcs(i,:)-vcs(j,:)
              dist1=SQRT(DOT_PRODUCT(vtp1,vtp1))
              IF(dist1<=plen(2) .AND.                                &
               & SQRT(DOT_PRODUCT(vcs(i,:),vcs(i,:)))<=plen(1) .AND. &
               & SQRT(DOT_PRODUCT(vcs(j,:),vcs(j,:)))<=plen(1)) THEN
                 np=np+1
                 pclst(np,1:3)=vcs(i,:)
                 pclst(np,4:6)=vcs(j,:)
                 pclst(np,7)=SQRT(DOT_PRODUCT(vcs(i,:),vcs(i,:)))
                 pclst(np,8)=SQRT(DOT_PRODUCT(vcs(j,:),vcs(j,:)))
                 pclst(np,9)=dist1
              END IF
           END DO
        END DO
        DO i=1,np
            CALL sort(pclst(i,7:8),2)
        END DO
        WRITE(o1,*) "----Pairs----"
        DO i=1,np
           WRITE(o1,110) (pclst(i,j),j=1,9)
           110 FORMAT(2X,2(3F8.4,1X),2F8.4,1X,F8.4)
        END DO
  !triplets
        nt=0
        DO i=1,nv-2
           DO j=i+1,nv-1
              DO k=j+1,nv
                 vtp1(:)=vcs(i,:)-vcs(j,:)
                 vtp2(:)=vcs(i,:)-vcs(k,:)
                 vtp3(:)=vcs(j,:)-vcs(k,:)
                 dist1=SQRT(DOT_PRODUCT(vtp1,vtp1))
                 dist2=SQRT(DOT_PRODUCT(vtp2,vtp2))
                 dist3=SQRT(DOT_PRODUCT(vtp3,vtp3))
                 IF(dist1<=tlen(2) .AND. dist2<=tlen(2) .AND.           &
                  & dist3<=tlen(2) .AND.                                &
                  & SQRT(DOT_PRODUCT(vcs(i,:),vcs(i,:)))<=tlen(1) .AND. &
                  & SQRT(DOT_PRODUCT(vcs(j,:),vcs(j,:)))<=tlen(1) .AND. &
                  & SQRT(DOT_PRODUCT(vcs(k,:),vcs(k,:)))<=tlen(1)) THEN
                    nt=nt+1
                    tclst(nt,1:3)=vcs(i,:)
                    tclst(nt,4:6)=vcs(j,:)
                    tclst(nt,7:9)=vcs(k,:)
                    tclst(nt,10)=SQRT(DOT_PRODUCT(vcs(i,:),vcs(i,:)))
                    tclst(nt,11)=SQRT(DOT_PRODUCT(vcs(j,:),vcs(j,:)))
                    tclst(nt,12)=SQRT(DOT_PRODUCT(vcs(k,:),vcs(k,:)))
                    tclst(nt,13)=dist1
                    tclst(nt,14)=dist2
                    tclst(nt,15)=dist3
                 END IF
              END DO
           END DO
        END DO
        DO i=1,nt
            CALL sort(tclst(i,10:12),3)
            CALL sort(tclst(i,13:15),3)
        END DO
        WRITE(o1,*) "----Triplets----"
        DO i=1,nt
           WRITE(o1,120) (tclst(i,j),j=1,15)
           120 FORMAT(2X,3(3F8.4,1X),3F8.4,1X,3F8.4)
        END DO
!
        END SUBROUTINE
!
        END MODULE clsts
