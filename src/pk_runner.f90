!=======================================================================================================================
!
!  *************
!  * pk_runner *
!  *************
!
!  Purpose:
!  --------
!  Holds subroutines required to pre-process, run, and post-procdess the point kinetics simulation.
!
!=======================================================================================================================
MODULE pk_runner
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       15/03/19   K. Luszczek              Original code
!
!=======================================================================================================================
!  Import subroutines:
   USE pk_methods     , ONLY: point_kinetics_step , pk_methods_init
   USE pk_time_control, ONLY: pk_time_step_size   , pk_time_save_check
   USE pk_read_input  , ONLY: rho_dollar          , process_input
   USE pk_write_output, ONLY: pk_write_to_file    , pk_plot_select
   USE pk_data_types  , ONLY: solution_storage_ini, solution_storage_add, solution_storage_process
!  Import variables:
   USE pk_kinds       , ONLY: dp
   USE pk_data_types  , ONLY: solution_at_step, gen_time, lambda_l, beta   , beta_l    , end_time, time_reset_max, &
                              min_delta_t     , set_data, a_table , b_table, time_table, pc_group, type_of_input
!
   IMPLICIT NONE
!
   PRIVATE
!
!  Public subroutines
!  ------------------
   PUBLIC :: run_point_kinetics ! Subroutine that performs one point kinetic time step.
   PUBLIC :: pre_processing     ! Subroutine to proces input files.
   PUBLIC :: post_processing    ! Subroutine to process the results of the run.
!
!  Module variables
!  ------------------
   TYPE(solution_at_step), SAVE :: bos, eos ! Information about solution at the beginning and at the end of a step.
   REAL(dp),               SAVE :: delta_t  ! Time-step size.
!
!=======================================================================================================================
CONTAINS
!=======================================================================================================================
!
!  ***********************
!  * run_point_kinetics *
!  ***********************
!
!
!=======================================================================================================================
SUBROUTINE run_point_kinetics()
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Local
!
   CHARACTER(120) :: reactivity_file_path  ! Reactivity input file path.
   CHARACTER(120) :: kin_param_file_path   ! Kinetic parameters input file path.
   INTEGER        :: repeat_counter = 0    ! Keeps track of how many times the time step was rejected.
   INTEGER        :: step_counter   = 0    ! Keeps track of how many time steps were calculated in total.
   INTEGER        :: ierror
   INTEGER        :: inner_iter = 0        ! Number of inner iterations inside point_kinetics_step().
   REAL(dp)       :: time_current = 0.0_dp ! Current total time in the simulation.
   REAL(dp)       :: delta_t_save          ! Temporarily hold the delta_t value prior to the new estimate. (s)
!   REAL(dp)       :: rho_dollar_eos        ! Reactivity at the end of step ($)
   LOGICAL        :: time_step_accept = .FALSE. ! Logical variable which indicates whether the time step is accepted.
!
!  Calculations
!
!  Main Loop
!  ---------
   main_simulation_loop: DO
!
      step_counter = step_counter + 1     ! Keep count of performed kinetic time-steps.
!
!  CALL the time step loop
!  -----------------------------
      CALL INNER_time_step_loop()
!
!  Update values for the next step
!  -----------------------------
      bos%c = eos%c
      bos%n2p = eos%n2p
      bos%n1p = eos%n1p
      eos%rho = eos%rho / beta ! Change the eos rho back to $ for later printing
!
   IF ( time_current >= end_time) EXIT ! Exit when the current time is larger than
!                                        the desired simulation time.
!
!  Check if save conditions are met.
!  If 'yes' then add the current step to the final solution.
!  -----------------------------
   CALL INNER_check_and_save()
!
   END DO main_simulation_loop
!
!
!  --------------------------------------------
!  Internal subroutines
!  --------------------------------------------
!
   CONTAINS
!=======================================================================================================================
   SUBROUTINE INNER_time_step_loop()
!=======================================================================================================================
      time_step_loop: DO repeat_counter = 1, time_reset_max
!
!        Determine the beginning of step rho_dollar
!        -----------------------------
         bos%rho = rho_dollar(&
            time_current  = time_current , & ! IN
            type_of_input = type_of_input, & ! IN
            a_table       = a_table      , & ! IN
            b_table       = b_table      , & ! IN
            time_table    = time_table   , & ! IN
            lambda        = lambda_l(1))   & ! IN, OPTIONAL
            * beta
!
            time_current = time_current + delta_t ! Increase the current simulation time by adding a time step
            IF ( time_current >= end_time) EXIT   ! Exit when the current time is larger than
!                                                   the desired simulation time.
!
!        Determine the end of step rho_dollar
!        ----------------------------- 
         eos%rho = rho_dollar(&
            time_current  = time_current , & ! IN
            type_of_input = type_of_input, & ! IN
            a_table       = a_table      , & ! IN
            b_table       = b_table      , & ! IN
            time_table    = time_table   , & ! IN
            lambda        = lambda_l(1))   & ! IN, OPTIONAL
            * beta
!
!        Calculate new n_eos and c_eos
!        -----------------------------
         CALL point_kinetics_step(&
            gen_time   = gen_time   , & ! IN
            lambda_l   = lambda_l   , & ! IN
            beta       = beta       , & ! IN
            beta_l     = beta_l     , & ! IN
            rho_eos    = eos%rho    , & ! IN
            rho_bos    = bos%rho    , & ! IN
            delta_t    = delta_t    , & ! IN
            n_bos      = bos%n2p    , & ! IN
            c_bos      = bos%c      , & ! IN
            n_eos      = eos%n2p    , & ! OUT
            n_eos_p1   = eos%n1p    , & ! OUT
            c_eos      = eos%c      , & ! OUT
            inner_iter = inner_iter)    ! OUT
!
!        Calculate the new time-step size (will be used in a reiteration or as the next step time-step).
!        Accept or reject the current time step.
!        -----------------------------
         delta_t_save = delta_t
         CALL pk_time_step_size(&
            delta_t          = delta_t          , & ! INOUT
            n_eos            = eos%n2p          , & ! IN
            n_eos_p1         = eos%n1p          , & ! IN
            time_current     = time_current     , & ! IN
            time_step_accept = time_step_accept)    ! OUT 
!
         IF ( time_step_accept ) THEN ! Check if the current time-step n_eos value was accepted
            EXIT
!        Break the loop if the step size is already at its minimum allowed value
         ELSE IF ( delta_t < min_delta_t ) THEN
            WRITE(*,*) 'WARNING ** Minimum allowed time-step reached at step: ', &
            step_counter, time_current
            EXIT
         ELSE
!           If the step was rejected do the following:
            IF (repeat_counter == time_reset_max ) THEN
!              If the maximum number of iterations was reached: display warning and do not reverse the time
               WRITE(*,*) ' WARNING ** Maximum number of time-step rejections reached at step: ', &
               step_counter, time_current
               EXIT
            ELSE
!              If the step will be re-iterated: reverse the time 
               time_current = time_current - delta_t_save
            ENDIF
         ENDIF
!
      END DO time_step_loop
!
   END SUBROUTINE INNER_time_step_loop
!=======================================================================================================================
   SUBROUTINE INNER_check_and_save()
!=======================================================================================================================

!     Locals
!
      LOGICAL :: save_check ! Logical flag to prompt the printout function.
!
!     Check if save conditions are met.
!     If 'yes' then add the current step to the final solution.
!     -----------------------------
      CALL pk_time_save_check(&
         time_current = time_current, & ! IN
         save_check   = save_check)     ! IN
!
      IF (save_check .OR. (repeat_counter > 1)) THEN
         CALL solution_storage_add(&
            IN_solution            = eos           , & ! IN
            IN_time_reject_count   = repeat_counter, & ! IN
            IN_implicit_iter_count = inner_iter    , & ! IN
            IN_step_number         = step_counter  , & ! IN
            IN_time_at_eos         = time_current )    ! IN
      END IF
!
   END SUBROUTINE
!
END SUBROUTINE run_point_kinetics
!=======================================================================================================================
!
!  ******************
!  * pre_processing *
!  ******************
!
!  Purpose: Pre-processing for point kientics program. Reads input files, inititlizes required modules.
!
!=======================================================================================================================
SUBROUTINE pre_processing()
!=======================================================================================================================
!
!  Locals
!
   CHARACTER(120) :: reactivity_file_path ! Reactivity input file path
   CHARACTER(120) :: kin_param_file_path  ! Kinetic parameters input file paths
   INTEGER  :: ierror
!
!  Actions
!
!  Read and process the input files
!  -------------------------------------------------------
   kin_param_file_path = 'reactor.in'
   reactivity_file_path = 'reactivity.in'
!
   CALL process_input(&
      kin_param_file_path  = kin_param_file_path , & ! IN
      reactivity_file_path = reactivity_file_path, & ! IN
      ierror               = ierror)                 ! OUT
!
!  Run initialization subroutines
!  -------------------------------------------------------
   CALL pk_methods_init()
   CALL pk_runner_init()
!
END SUBROUTINE pre_processing
!=======================================================================================================================
!
!  *******************
!  * post_processing *
!  *******************
!
!  Purpose: Post-processing of the results. Writing out to file(s), creating selected plots.
!
!=======================================================================================================================
SUBROUTINE post_processing()
!=======================================================================================================================
!
!  Actions
!
   CALL solution_storage_process()
   CALL pk_write_to_file()
   CALL pk_plot_select(plot='PLOT_N_RHO')
   CALL pk_plot_select(plot='PLOT_C')
!
END SUBROUTINE post_processing
!=======================================================================================================================
!
!  ******************
!  * pk_runner_init *
!  ******************
!
!  Purpose: Subroutine initiating the point kinetics runner.
!
!=======================================================================================================================
SUBROUTINE pk_runner_init()
!=======================================================================================================================
! Record of revisions:
!        Date       Programmer               Description of changes
!        -----------------------------------------------------------
!        27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Locals
!
   INTEGER:: idx ! Loop index
!
!  Allocate memory for precursor density arrays
!  -------------------------------------------------------
   ALLOCATE(bos%c(pc_group))
   ALLOCATE(eos%c(pc_group))
!
!  Calculate the initial precursor densities
!  -------------------------------------------------------
   initialc_bos: DO idx = 1, pc_group
      bos%c(idx) = beta_l(idx) / ( lambda_l(idx) * gen_time)
   ENDDO initialc_bos
   eos%c = bos%c
!
!  Initalize the linked list (for storing the solution at each calculated step) by passing the first entry.
!  -------------------------------------------------------
   CALL solution_storage_ini(&
      IN_solution            = bos   , & ! IN
      IN_time_reject_count   = 0     , & ! IN
      IN_implicit_iter_count = 0     , & ! IN
      IN_step_number         = 0     , & ! IN
      IN_time_at_eos         = 0.0_dp )  ! IN
!
!  Set the initial time step to the minimum allowed value
!  -------------------------------------------------------
   delta_t = min_delta_t
!
END SUBROUTINE pk_runner_init
!=======================================================================================================================
END MODULE pk_runner
!=======================================================================================================================
