!=======================================================================================================================
!
!  *******************
!  * pk_time_control *
!  *******************
!
!  Purpose:
!  --------
!  This module contains routines needed for point kinetic time-step size control.
!
!  Usage:
!  ------
!
!  CALL pk_time_step_size   to perform the time-step size adjustment.
!  CALL pk_time_save_check  to check if the current time step should be saved for post-processing
!
!=======================================================================================================================
MODULE pk_time_control
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       07/01/19   K. Luszczek              Original code
!
!=======================================================================================================================
!  Import subroutines
   USE pk_write_output , ONLY: pk_write_to_file
   USE pk_data_types   , ONLY: set_data, make_knot_table
!  Import variables
   USE pk_kinds        , ONLY: dp
   USE pk_write_output , ONLY: pk_write_to_file
   USE pk_data_types   , ONLY: min_ratio, max_ratio, tol_rel, edit_delta_t, min_delta_t , knot_table, end_time
!
   IMPLICIT NONE
!
   PRIVATE
!
   REAL(dp), SAVE :: last_recom_dt    = 0.0_dp  ! Saves last recommended time-step size. May be different than the 
!                                                 one actually used in the calculations.
   REAL(dp), SAVE :: next_time_to_print         ! Next time step that should be printed out. Can't be omitted
   INTEGER,  SAVE :: knot_table_index = 1       ! Keeps track of the last used entry from the knot_table_index
   REAL(dp), SAVE :: delta_t_avg      = 0.0_dp  ! An average time step
!
!  Public subroutines and variables
!  --------------------------------
   PUBLIC :: pk_time_control_init ! Initializes the module
   PUBLIC :: pk_time_step_size    ! Subroutine that performs the time-step size adjustment
   PUBLIC :: pk_time_save_check   ! Check if the results should be saved to the linked list.
!
!=======================================================================================================================
CONTAINS
!=======================================================================================================================
!
!  *********************
!  * pk_time_step_size *
!  *********************
!
!  Purpose: Return the recommended time step.
! 
!=======================================================================================================================
SUBROUTINE  pk_time_step_size(&
   delta_t,                   & ! INOUT (s)
   n_eos,                     & ! IN
   n_eos_p1,                  & ! IN
   time_current,              & ! IN
   time_step_accept)            ! IN
!=======================================================================================================================
! Record of revisions:
!       Date     Programmer  Description of change
!       ====     ==========  =====================
!       08/01/19 K. Luszczek Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(INOUT) :: delta_t          ! Kinetic time-step. Passed on as current, returned as new (s)
   REAL(dp), INTENT(IN)    :: n_eos            ! End of step neutron density   
   REAL(dp), INTENT(IN)    :: n_eos_p1         ! End of step neutron density from the second order method
   REAL(dp), INTENT(INOUT) :: time_current     ! Current time-stamp (s)
   LOGICAL,  INTENT(OUT)   :: time_step_accept ! To return true or false back to the calling program
!
!  Locals
!
   REAL(dp) :: lte    ! Local truncation error estimate
   REAL(dp) :: tol    ! Requested tolerance
   REAL(dp) :: dt_rec ! Recommended new time step size
   REAL(dp) :: safety = 0.9_dp ! Safety factor
!
   lte = ABS(n_eos - n_eos_p1) ! Estimate the local trun. error
   tol = n_eos * tol_rel
!
!  Check whether the time step is to be accepted or rejected
!  --------------------------------
   CALL INNER_acceptance_check()
!
   CALL INNER_recommend_step()
!
   CALL INNER_adjust_step_for_edit()
!
   delta_t = dt_rec
!
!  --------------------------------------------
!  Internal subroutines
!  --------------------------------------------
!
CONTAINS
!=======================================================================================================================
   SUBROUTINE INNER_acceptance_check()
!=======================================================================================================================
      IF ( lte <= tol ) THEN
         time_step_accept = .TRUE.
      ELSE
         time_step_accept = .FALSE.
      ENDIF
   END SUBROUTINE INNER_acceptance_check
!=======================================================================================================================
   SUBROUTINE INNER_recommend_step()
!=======================================================================================================================
!     LTE might be very close to zero if n_eos -> n_bos. In that case set to some small non-zero value to avoid numerical problems.
      IF ( lte <= EPSILON(lte)*10 ) THEN 
         lte = EPSILON(lte)*100     
         WRITE(*,*) 'min lte !! '    
      ENDIF
!     Recommend the new time step 
!     The new time step HAS TO BE smaller than 1.0! 
      dt_rec = MAX(delta_t * MIN(MAX(min_ratio,SQRT(tol/lte)*safety),max_ratio), min_delta_t)
   END SUBROUTINE INNER_recommend_step
!=======================================================================================================================
   SUBROUTINE INNER_adjust_step_for_edit()
!=======================================================================================================================
!  Check if the edit time step will not be hopped over in the next time_step_size. Add 25% to the current time step.
   IF ((time_current+1.25_dp*dt_rec)>(next_time_to_print)) THEN
!     Check if the next edit time and current time are identical, or close enough that they are indistiguishable
!     with the computer's epsilon. (more or less a proper way of comparing float numbers??)
!     If yes, set the current time to be identical to the next edit time. This avoids some numerical issues down the line.
!     (Previously it was inadvertently bringing the time step to the minimum allowed value)
      IF (ABS(next_time_to_print - time_current)/next_time_to_print < epsilon(time_current)*10) THEN
         time_current = next_time_to_print
!     With that change, check if the recommended time step is bigger than edit_delta_t. If it is you will accidentally jump
!     over the step. In that case, limit the size. This may artificially limit the time step size.
         IF ((knot_table_index+1) > SIZE(knot_table)) THEN
            dt_rec = min_delta_t
         ELSE IF (dt_rec >= (knot_table(knot_table_index+1)-knot_table(knot_table_index))) THEN
            dt_rec = knot_table(knot_table_index+1)-knot_table(knot_table_index)
         ENDIF
      ELSE 
!        If the distance from the current time point to the next edit time point is larger than
!        the computer's espilon:a
!        -adjust the time step to fall exactly on the edit time step 
         dt_rec = next_time_to_print - time_current
!        Limit the smallest dt_rec size to be at least the minimum time step size
         IF (dt_rec < min_delta_t) THEN
            dt_rec = min_delta_t
         ENDIF
      END IF
   END IF
   END SUBROUTINE INNER_adjust_step_for_edit
END SUBROUTINE pk_time_step_size
!=======================================================================================================================
!
!  **********************
!  * pk_time_save_check *
!  **********************
!
!  Purpose: Return the recommended time step.
! 
!=======================================================================================================================
SUBROUTINE  pk_time_save_check(&
   time_current,               & ! IN
   save_check)                   ! IN
!=======================================================================================================================
! Record of revisions:
!       Date     Programmer  Description of change
!       ====     ==========  =====================
!       08/01/19 K. Luszczek Original code 
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: time_current ! Current time-step (s)
   LOGICAL,  INTENT(OUT):: save_check   ! Flag to prompt printout function
!
!   Calculations
!
   IF (time_current >= next_time_to_print .OR. &
    abs(time_current-next_time_to_print)/time_current < EPSILON(time_current)*10) THEN
!
      save_check = .TRUE.
      knot_table_index = knot_table_index + 1
      IF (knot_table_index > SIZE(knot_table)) THEN
         next_time_to_print = end_time
      ELSE 
         next_time_to_print = knot_table(knot_table_index)
      END IF
   ELSE
      save_check = .FALSE.
   END IF 
!  
END SUBROUTINE pk_time_save_check
!=======================================================================================================================
!
!  ************************
!  * pk_time_control_init *
!  ************************
!
!  Purpose: Initializes the module
! 
!=======================================================================================================================
SUBROUTINE  pk_time_control_init()
!=======================================================================================================================
! Record of revisions:
!       Date     Programmer  Description of change
!       ====     ==========  =====================
!       08/01/19 K. Luszczek Original code 
!
!=======================================================================================================================
!
   CALL make_knot_table()
   next_time_to_print = knot_table(knot_table_index) ! Initialize the as the first entry in the knot table
   knot_table_index = knot_table_index
! 
END SUBROUTINE pk_time_control_init
!=======================================================================================================================
END MODULE pk_time_control
!=======================================================================================================================