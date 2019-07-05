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
   USE pk_data_types   , ONLY: set_data
!  Import variables
   USE pk_kinds        , ONLY: dp
   USE pk_write_output , ONLY: pk_write_to_file
   USE pk_data_types   , ONLY: min_ratio, max_ratio, eps_rel, edit_delta_t, min_delta_t
!
   IMPLICIT NONE
!
   PRIVATE
!
   REAL(dp), SAVE :: last_recom_dt      = 0.0_dp  ! Saves last recommended time-step size. May be different than the 
!                                                   one actually used in the calculations.
   REAL(dp), SAVE :: next_time_to_print = 0.0_dp  ! Next time step that should be printed out. Can't be omitted
   REAL(dp), SAVE :: delta_t_avg        = 0.0_dp  ! An average time step
!
!  Public subroutines and variables
!  --------------------------------
   PUBLIC :: pk_time_step_size  ! Subroutine that performs the time-step size adjustment
   PUBLIC :: pk_time_save_check ! Check if the results should be saved to the linked list.
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
   REAL(dp), INTENT(INOUT) :: delta_t        ! Kinetic time-step. Passed on as current, returned as new (s)
   REAL(dp), INTENT(IN) :: n_eos             ! End of step neutron density   
   REAL(dp), INTENT(IN) :: n_eos_p1          ! End of step neutron density from the second order method
   REAL(dp), INTENT(IN) :: time_current      ! Current time-stamp (s)
   LOGICAL,  INTENT(OUT) :: time_step_accept ! To return true or false back to the calling program
!
!  Locals
!
   REAL(dp) :: lte    ! Local truncation error estimate
   REAL(dp) :: tol    ! Requested tolerance
   REAL(dp) :: dt_rec ! Recommended new time step size
   REAL(dp) :: safety = 0.9_dp ! Safety factor
!
   lte = ABS(n_eos - n_eos_p1) ! Estimate the local trun. error
   tol = n_eos * eps_rel
!
!  Check whether the time step is to be accepted or rejected
!
   IF ( lte <= tol ) THEN
      time_step_accept = .TRUE.
   ELSE
      time_step_accept = .FALSE.
   ENDIF
!
   IF ( lte <= 1.0e-14_dp ) THEN  ! LTE might be very close to zero if n_eos -> n_bos. In that case set to some small
      lte = 1.0e-14_dp            ! non-zero value to avoid numerical problems.
   ENDIF
!
   dt_rec = delta_t * MIN(MAX(min_ratio,SQRT(tol/LTE)*safety),max_ratio)    ! Recommend the new time step
!
   IF (dt_rec < min_delta_t) THEN  ! Check if smaller than the minimum allowed time step
      dt_rec = min_delta_t 
   END IF
!
! Check if the edit time step will not be hopped over in the next time_step_size
   IF ((time_current+dt_rec)>next_time_to_print) THEN
      IF ((next_time_to_print - time_current) < min_delta_t) THEN
         ! The distance from the current time point to the next edit time point is smaller
         ! than the minimum allowable time step. Set the time step to minimum.
         ! Set the recommended time step to the minimum allowable step.
         dt_rec = min_delta_t
      ELSE 
         ! If the distance from the current time point to the next edit time point is larger than
         ! the minimum allowable time step. Adjust the time step to fall exactly on the 
         ! edit time step 
         dt_rec = next_time_to_print - time_current
      END IF
   END IF
!
delta_t = dt_rec
!
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
   IF (time_current >= next_time_to_print) THEN
      save_check = .TRUE.
      next_time_to_print = next_time_to_print + edit_delta_t
!
   ELSE
      save_check = .FALSE.
!
   END IF 
!  
END SUBROUTINE pk_time_save_check
!=======================================================================================================================
END MODULE pk_time_control
!=======================================================================================================================
