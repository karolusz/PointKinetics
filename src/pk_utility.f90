!=======================================================================================================================
!
!  **************
!  * pk_utility *
!  **************
!
!  Purpose:
!  --------
!  Module containing some utility functions.
!
!  Usage:
!  ------
!  CALL pk_sort  to sort 1D array (limited to an array of _dp kind)
!
!=======================================================================================================================
MODULE pk_utility
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       12/08/19   K. Luszczek              Original code.
!
!=======================================================================================================================
!   Import variables
    USE pk_kinds, ONLY: dp
!
    IMPLICIT NONE
!
    PRIVATE
!
!   Public subroutines
!   ------------------
    PUBLIC :: pk_hpsort_eps ! Subroutine which process the input files.
    PUBLIC :: pk_remove_duplicates ! Removes duplicate entires from the table.
!
!=======================================================================================================================
CONTAINS
!=======================================================================================================================
SUBROUTINE pk_hpsort_eps(&    
    array_to_sort ) ! IN
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       12/08/19   K. Luszczek              Original code.
!
!=======================================================================================================================
!
!  Arguments
!
    REAL(dp), INTENT(INOUT) :: array_to_sort(:)
!
!   Locals
!
    INTEGER  :: n   ! Size of the array to be sorted
    REAL(dp) :: eps ! Machine espislon for selected type
!
    n = SIZE(array_to_sort)
    eps = EPSILON(array_to_sort(1))
!
    CALL hpsort_eps(&
        n = n,             & ! IN
        eps = eps,         & ! IN
        ra = array_to_sort)  ! INOUT
END SUBROUTINE pk_hpsort_eps
!=======================================================================================================================
SUBROUTINE hpsort_eps(&
    n,   & ! IN 
    eps, & ! INOUT
    ra)    ! IN
!=======================================================================================================================
!
!   Arguments
!
    INTEGER , INTENT(IN)    :: n     ! Size of the array to be sorted
    REAL(dp), INTENT(IN)    :: eps   ! Machine epsilon for dp kind
    REAL(dp), INTENT(INOUT) :: ra(n) ! Array to be sorted
!
!   Locals
!
    INTEGER :: i
    INTEGER :: ir
    INTEGER :: j
    INTEGER :: l
    INTEGER :: iind
    INTEGER :: ind(n) ! Index array
    REAL(dp) :: rra
!
    DO i =1, n
        ind(i) = i
    END DO
!
  ! nothing to order
    IF (n.lt.2) return  
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1  
  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 
!
!  internal function 
!  compare two real number and return the result
!
    LOGICAL FUNCTION hslt( a, b )
    REAL(DP) :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
    END FUNCTION hslt
!
END SUBROUTINE hpsort_eps
!=======================================================================================================================
SUBROUTINE pk_remove_duplicates(&    
    array_to_sort ) ! IN
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       12/08/19   K. Luszczek              Original code.
!
!=======================================================================================================================
!
!  Arguments
!
    REAL(dp), ALLOCATABLE, INTENT(INOUT) :: array_to_sort(:)
!
!   Locals
!
    INTEGER               :: n   ! Size of the array to be sorted
    INTEGER               :: idx ! Index variable
    INTEGER               :: idy ! Index variable
    INTEGER, ALLOCATABLE  :: dup(:)         ! Array of indecies with duplicates
    INTEGER               :: dup_size
    REAL(dp)              :: eps            ! Machine espislon for selected type
    REAL(dp), ALLOCATABLE :: no_dups_arr(:) ! Array without duplicates
!
    n = SIZE(array_to_sort)
    eps = EPSILON(array_to_sort(1))
!
    ALLOCATE(dup(n))
    dup(1) = 1
    DO idx = 1, n-1
        IF(abs(array_to_sort(idx)-array_to_sort(idx+1)) <= eps) THEN
            dup(idx+1) = 0
        ELSE 
            dup(idx+1) = 1
        END IF
    END DO
!   
    dup_size = SUM(dup)
    ALLOCATE(no_dups_arr(dup_size))
!
    idy = 1
    DO idx = 1, n
        IF (dup(idx) > 0) THEN
            no_dups_arr(idy) = array_to_sort(idx)
            idy = idy + 1
        ENDIF
    END DO
!
    CALL MOVE_ALLOC(no_dups_arr,array_to_sort)
!
END SUBROUTINE pk_remove_duplicates
!=======================================================================================================================
END MODULE pk_utility
!=======================================================================================================================