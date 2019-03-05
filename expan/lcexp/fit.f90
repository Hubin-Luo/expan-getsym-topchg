        SUBROUTINE fit()
!----------------------------------------------
!A least-square fit                    HBL
!----------------------------------------------
        USE ctrl; USE cellpos; USE abin
        IMPLICIT NONE
        INTEGER :: i,j,k,m,icf,iter
        REAL(q) :: rme,me,smsq
        REAL(q) :: fen(ne)
        CHARACTER(3) :: fmt 
        REAL(q),ALLOCATABLE :: fun(:,:),aa(:,:)
        REAL(q),ALLOCATABLE :: eci(:),enf(:)     !dot product of eng.fun
!
!allocate eci and enf
!
        ALLOCATE(eci(nterm))
        ALLOCATE(enf(nterm))
!
!allocate the matrices
!
        ALLOCATE(fun(ne,nterm))
        ALLOCATE(aa(nterm,nterm))
!
!assign and calculate the matrices
!
        fun(1:ne,1)=1._q
        IF(ncomp==2) THEN
           DO i=1,ne
              DO j=2,nterm
                 IF(j<=is+1) THEN
                    fun(i,j)=pos(ist(i))%fg1(j-1,1)
                 ELSE IF(j<=is+ip+1) THEN
                    fun(i,j)=pos(ist(i))%fg2(j-1-is,1)
                 ELSE
                    fun(i,j)=pos(ist(i))%fg3(j-1-is-ip,1)
                 END IF
              END DO
           END DO
        END IF
!
!        IF(ncomp==3) THEN
!           DO i=1,ne
!              DO j=2,nterm
!                 IF(j<=2*is+1) THEN
!                    k=INT( (j-1+1.5_q)/2 )
!                    m=j-2*INT((j-1-0.5_q)/2)
!                    fun(i,j)=pos(ist(i))%fg1(k,m)
!                 ELSE IF(j<=2*is+4*ip+1) THEN
!                    k=INT( (j-1-2*is+1.5_q)/4 )
!                    m=j-1-2*is-4*INT( (j-1-2*is-0.5_q)/4 )
!                    fun(i,j)=pos(ist(i))%fg2(k,m)
!                 ELSE
!                    k=INT( (j-1-2*is-4*ip+1.5_q)/8 )
!                    m=j-1-2*is-4*ip-8*INT( (j-1-2*is-4*ip)/8 )
!                    fun(i,j)=pos(ist(i))%fg3(k,m)
!                 END IF
!              END DO
!           END DO
!        END IF
!
        IF(ncomp==3) THEN
          DO i=1,ne
             icf=1
             DO k=1,is         !transfer the single point cluster functions
                DO m=1,2
                   icf=icf+1
                   fun(i,icf)=pos(ist(i))%fg1(k,m)
                END DO
             END DO
             DO k=1,ip         !transfer the pair cluster functions
                DO m=1,4
                   icf=icf+1
                   fun(i,icf)=pos(ist(i))%fg2(k,m)
                END DO
             END DO
             DO k=1,it         !transfer the triplet cluster functions
                DO m=1,8
                   icf=icf+1
                   fun(i,icf)=pos(ist(i))%fg3(k,m)
                END DO
             END DO
          END DO
        END IF
!
        FORALL(i=1:nterm,j=1:nterm)
           aa(i,j)=DOT_PRODUCT(fun(:,i),fun(:,j))
        END FORALL
        DO i=1,nterm
           enf(i)=DOT_PRODUCT(eng(:),fun(:,i))
        END DO
        CALL sor(aa,enf,nterm,1.2_q,iter,eci)
        WRITE(o3,*) "*************************"
        WRITE(o3,*) "Information from fitting:"
        WRITE(o3,*) "*************************"
        WRITE(o3,"(/,1X,'Stopped at: ',I5)") iter
        fmt="???"
        WRITE(fmt,'(I3)') nterm
        WRITE(o3,"(1X,'ECIs: ',"//fmt//"(E16.8)/)") &
              &(eci(i),i=1,nterm)
!
! calculate the root mean square deviation
!
        fen=0._q; smsq=0._q; rme=0._q; me=0._q
        DO i=1,ne
           DO j=1,nterm
              fen(i)=fen(i)+eci(j)*fun(i,j)
           END DO
           rme=rme+(fen(i)-eng(i))*(fen(i)-eng(i))
        END DO
        rme=SQRT(rme/ne)
        DO i=1,ne
           me=me+ABS(fen(i)-eng(i))
        END DO
        me=me/ne
        WRITE(o3,"(1X,'Standard deviation: ', E10.2)") rme
        WRITE(o3,"(1X,'Mean deviation: ', E10.2)") me
        WRITE(o3,*) "Input and fitted energies:"
        WRITE(o3,"(1X,'Esub= ',I8)") nesub
        DO i=1,ne
           WRITE(o3,"(1X,I4,2(1X,F15.6))") ist(i),eng(i),fen(i)
        END DO
!
        DEALLOCATE(fun,aa,enf,eci); DEALLOCATE(ist,eng)
!
        END SUBROUTINE


