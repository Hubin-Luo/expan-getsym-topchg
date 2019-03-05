!**********************************************************************
!This subroutine uses the spglib to get the symmetry operation matrices
!The output arrays "mrot" and "vtran" store the rotation matrices and
!tranlation vectors. They are not allocated before provoking but will
!be allocated after. You should input the precision of finding symmetry
!which is usually 1.e-5
! BE CAREFUL:
! The number of operations is controlled by the last dimension. The matrix 
! elements are stored in column order. So, if you want to do symmetry
! operation to a vector, 
! use MATMUL(vec, mrot(:,:,iop)) + vtran(:,iop) to get the new one.
!
! Compilation needs the interface in spglib_f08.f90 to link to 
! libsymspg.a
!---------------------------------------------------------------H. B. Luo
!
        MODULE symmat
!
        INTEGER :: ntyp, natom, nop
        INTEGER, ALLOCATABLE :: ispec(:)
        INTEGER, ALLOCATABLE :: mrot(:,:,:)
        REAL(8), ALLOCATABLE :: vtran(:,:), pos(:,:)
        REAL(8) :: latt(3,3)
!
        CONTAINS
!
        SUBROUTINE getsym(symprec)
!        
        USE spglib_f08
        IMPLICIT NONE
!
        INTEGER :: i, j, k, Ispg, max_size
        INTEGER :: iStrt, iEnd
        INTEGER :: itmp(3,3)
        INTEGER, ALLOCATABLE :: rot(:,:,:)
        REAL(8), ALLOCATABLE :: trans(:,:)
        INTEGER, AlLOCATABLE :: itype(:)
        REAL(8) :: scalor, symprec
!
        max_size=5000
        ALLOCATE(rot(3,3,max_size)); ALLOCATE(trans(3,max_size))
        rot=0
!
        OPEN(unit=8, file='POSCAR', status='old')
        READ(8,*) ntyp
!
        READ(8,*) scalor
        DO j=1,3
           READ(8,*) (latt(i,j),i=1,3)
        END DO
        latt=scalor*latt
!
        ALLOCATE(ispec(ntyp))
        READ(8,*)
        READ(8,*) (ispec(i),i=1,ntyp)
!
        natom=0
        DO i=1,ntyp
           natom=natom+ispec(i)
        END DO
!
        ALLOCATE(itype(natom))
        iStrt=1; iEnd=0
        DO i=1, ntyp
           iEnd = iStrt + ispec(i) - 1
           itype(iStrt:iEnd) = i
           iStrt = iEnd + 1
        END DO
!
        ALLOCATE(pos(3,natom)) 
        READ(8,*)
        DO j=1, natom
           READ(8,*) (pos(i,j), i=1,3)
        END DO
!       WRITE(*,*) "Positions:"
!       DO j=1, natom
!          WRITE(*,*) (pos(i,j), i=1,3)
!       END DO
        CLOSE(8)
!
        Ispg=spg_get_symmetry(rot, trans, max_size, &
                &  latt, pos, itype, natom, symprec)
        IF(Ispg==0) THEN
          WRITE(*,*) "Failed to get symmetry"
        ELSE
          DO i=1, max_size
             itmp=rot(:,:,i)
             IF (MAXVAL(itmp*itmp)==0) THEN
                nop = i - 1
                EXIT
             END IF
          END DO
          ALLOCATE(mrot(3,3,nop)); ALLOCATE(vtran(3,nop))
          mrot=rot(:,:,1:nop)
          vtran=trans(:,1:nop)
          OPEN(9,file='SYMMAT',status='replace')
          DO i=1, nop
             WRITE(9,*) "Rotation Matrix"
             DO j=1,3
               WRITE(9,'(1X,3I4)') (mrot(k,j,i), k=1,3)
             END DO
             WRITE(9,*) "translation"
             WRITE(9,'(1X,3F10.6)') (vtran(k,i), k=1,3)
             WRITE(9,*)
          END DO
          CLOSE(9)
        END IF 
        DEALLOCATE(rot); DEALLOCATE(trans)
        DEALLOCATE(itype)
!
        END SUBROUTINE
!
        END MODULE
