!
        MODULE ctrl
          INTEGER,PARAMETER:: q=SELECTED_REAL_KIND(10)
          INTEGER,PARAMETER:: i1=7, o1=8, o2=9,o3=10
          REAL(q), PARAMETER:: thd=0.000001
          REAL(q), PARAMETER:: hhd=0.00000001
          CHARACTER(LEN=20) :: fnm
          INTEGER :: ncomp,nterm
        END MODULE
!
        MODULE abin
          USE ctrl
          INTEGER :: ne,nesub
          INTEGER,ALLOCATABLE :: ist(:)
          REAL(q),ALLOCATABLE :: eng(:)
        END MODULE
