        SUBROUTINE rotardaxis(axis, tranv, thetadeg, vecin, vecout, pov)
!---------------------------------------------------------------
! Rotate a vector around a given axis
! tranv is the vector to translate the axis to origin(not unique)
! The last parameter defines rotation of position(0) or vector(1)
!                                                      H. B. Luo
!---------------------------------------------------------------
        IMPLICIT NONE
        INTEGER :: i, j, k
        INTEGER,INTENT(IN) :: pov
        REAL(8),INTENT(IN) :: thetadeg
        REAL(8),INTENT(IN) :: axis(3), tranv(3), vecin(3)
        REAL(8),INTENT(OUT) :: vecout(3)
        REAL(8) :: vectmpin(4), vectmpout(4)
        REAL(8) :: rotxp(4,4), rotyp(4,4), rotzp(4,4), tranmp(4,4)
        REAL(8) :: rotxn(4,4), rotyn(4,4), rotfin(4,4),tranmn(4,4)
        REAL(8) :: sinalph, cosalph, sinbeta, cosbeta, theta, pi
!
        pi=ACOS(-1.d0)
!        axis=(/1.d0,1.d0,0.d0/); theta=pi
!        vecin=(/1.d0, 0.d0, 0.d0/)
        rotxp=0.d0; rotyp=0.d0; rotzp=0.d0; tranmp=0.d0; rotfin=0.d0
!
! set the rotation sin and cos of the axis around x and y. remember if the
! axis is x axis, we have to set sinalph and cosalph manually to avoid singularity
!
        IF(ABS(axis(2))<1.d-6 .AND. ABS(axis(3))<1.d-6) THEN
           cosalph=1.d0
           sinalph=0.d0
        ELSE
           sinalph=axis(2)/(SQRT(axis(2)*axis(2)+axis(3)*axis(3)))
           cosalph=axis(3)/(SQRT(axis(2)*axis(2)+axis(3)*axis(3)))
        END IF
        sinbeta=axis(1)/SQRT(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3))
        cosbeta=SQRT(axis(2)*axis(2)+axis(3)*axis(3))/ &
       &        SQRT(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3))
        theta=thetadeg/180*pi
!
        rotxp(1,1)=1.d0;
        rotxp(2,2)=cosalph
        rotxp(2,3)=-sinalph
        rotxp(3,2)=-rotxp(2,3)
        rotxp(3,3)=rotxp(2,2)
        rotxp(4,4)=1.d0
        rotxn=rotxp; rotxn(2,3)=-rotxn(2,3); rotxn(3,2)=-rotxn(3,2)
!
        rotyp(1,1)=cosbeta
        rotyp(1,3)=sinbeta
        rotyp(3,1)=-rotyp(1,3)
        rotyp(3,3)=rotyp(1,1)
        rotyp(2,2)=1.d0
        rotyp(4,4)=1.d0
        rotyn=rotyp; rotyn(1,3)=-rotyn(1,3); rotyn(3,1)=-rotyn(3,1)
!
        rotzp(1,1)=COS(theta)
        rotzp(1,2)=-SIN(theta)
        rotzp(2,1)=-rotzp(1,2)
        rotzp(2,2)=rotzp(1,1)
        rotzp(3,3)=1.d0
        rotzp(4,4)=1.d0
!
        DO i=1, 4
           tranmp(i, i)=1.d0
           tranmn(i, i)=1.d0
        END DO
        tranmn(1:3, 4)=-tranv
        tranmp(1:3, 4)=tranv
!
        vectmpin(1:3)=vecin
        IF (pov==0) THEN
            vectmpin(4)=1.d0     ! rotate position
        ELSE IF (pov==1) THEN
            vectmpin(4)=0.d0     ! rotate vector
        ELSE
            WRITE(*,*) "[ROTARDAXIS]: Please define the rotation, position or vector"
            STOP
        END IF
        DO i=1, 4
           rotfin(i, i)=1.d0
        END DO
        rotfin=MATMUL(rotfin, tranmp)
        rotfin=MATMUL(rotfin, rotxn)
        rotfin=MATMUL(rotfin, rotyp)
        rotfin=MATMUL(rotfin, rotzp)
        rotfin=MATMUL(rotfin, rotyn)
        rotfin=MATMUL(rotfin, rotxp)
        rotfin=MATMUL(rotfin, tranmn)
!
        vectmpout=MATMUL(rotfin, vectmpin)
        vecout=vectmpout(1:3)
!
        END SUBROUTINE
