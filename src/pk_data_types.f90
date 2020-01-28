!=======================================================================================================================
!
!  *****************
!  * pk_data_types *
!  *****************
!
!  Purpose:
!  --------
!
!  Defines the derived data types and default parameters' values needed for the point kinetics module.
!
!  Usage:
!  ------
!  CALL data_set()              to change a value of a protected variable in this module.
!  CALL allocate_1d()           to allocate memory for the protected array in this module.
!  CALL deallocate_1d()         to deallocate the protected array in this module.
!  CALL solution_storage_ini()  to initialize the linked list of point kinetics results.
!  CALL solution_storage_add()  to add an entry to the solution_storage linked list.
!  CALL step_size_store_ini()   to initialize the linked list for storing time step sizes.
!  CALL step_size_store_add()   to add an entry to the step_size_store linked list.
!  CALL step_size_fetch()       to get the next entry in the step_size store linked list.
!
!=======================================================================================================================
MODULE pk_data_types
!=======================================================================================================================
! Record of revisions:
!       Date     Programmer  Description of change
!       ====     ==========  =====================
!       18/03/19 K. Luszczek Original code
!
!=======================================================================================================================
   USE pk_kinds      , ONLY: dp
   USE pk_enumeration, ONLY: C_FULLY_IMPLICIT    , C_CONSTANT_APPROX, C_LINEAR_APPROX , C_EXPONENT_APPROX, &
                             N_NO_FREQUENCY_TRANS, N_FREQUENCY_TRANS, DT_ADAPT        , DT_CONST         , &
                             OMEGA_YES           , OMEGA_NO         , CONV_N          , CONV_W
   USE pk_utility    , ONLY: pk_hpsort_eps       , pk_remove_duplicates
!
   IMPLICIT NONE
!
   PRIVATE
!
!  Public types
!  ------------------
   PUBLIC :: solution_at_step, solution_storage
!
!  Public variables
!  ------------------
   PUBLIC :: pc_group      , gen_time      , beta          , lambda_l      , beta_l
   PUBLIC :: n_solution    , pre_approx    , epsilon_n_eos , max_n_eos_iter, omega_option, tol_rel  , converg_opt
   PUBLIC :: end_time      , time_reset_max, edit_table    , edit_delta_t  , min_ratio   , max_ratio, min_delta_t
   PUBLIC :: dt_option     , dt_user       , knot_table
   PUBLIC :: a_table       , b_table       , time_table_rho, type_of_input
   PUBLIC :: sol_list_head , sol_list_tail
   PUBLIC :: fs_time_stamps, fs_n          , fs_rho        , fs_c
   PUBLIC :: exec_time  
!
!  Public subroutines
!  ------------------
   PUBLIC :: set_data
   PUBLIC :: set_dt_options
   PUBLIC :: allocate_1d, deallocate_1d
   PUBLIC :: solution_storage_ini, solution_storage_add, solution_storage_process
   PUBLIC :: step_size_store_ini , step_size_store_add , step_size_fetch
   PUBLIC :: make_knot_table
!
!  Interface for set functions
!  ------------------
   INTERFACE set_data
      MODULE PROCEDURE set_data_real_single
      MODULE PROCEDURE set_data_real_array
      MODULE PROCEDURE set_data_int
      MODULE PROCEDURE set_data_chr
   END INTERFACE
!
!  Derived data types
!  ------------------
TYPE :: solution_at_step
   REAL(dp) :: n2p = 1.0_dp ! Neutron (relative) density. Second-order method. (no unit). Total flux.
   REAL(dp) :: n1p = 1.0_dp ! Neutron (relative) density. First-order method. (no unit). Total flux
   REAL(dp) :: n2p_tsc = 1.0_dp ! Second-order neturon density or freq. trans. density for time control.
   REAL(dp) :: n1p_tsc = 1.0_dp ! First-order neutron denisty or freq. trans. denisty for time control.
   REAL(dp) :: rho = 0.0_dp ! Reactivity at a given step. ($)
   REAL(dp), ALLOCATABLE :: c(:) ! Delayed neutron precursor concentration. (no unit)
END TYPE solution_at_step
!
TYPE :: solution_storage
   TYPE(solution_storage), POINTER :: p_next_entry  ! Points to the subsequent entry in the linked list.
   TYPE(solution_at_step)          :: kinetics_data ! Stores point kinetic solution at end of step.
   INTEGER  :: time_reject_count   = 0      ! Number of times the step size was rejected during step calculation.
   INTEGER  :: implicit_iter_count = 0      ! Implicit iterations performed.
   INTEGER  :: step_number         = 0      ! The point-kinetic step number.
   REAL(dp) :: time_at_eos         = 0.0_dp ! The time at the end of step. (s)
END TYPE solution_storage
!
TYPE :: step_size_store
   TYPE(step_size_store), POINTER :: p_next_entry ! Points to the subsequent entry in the linked list.
   INTEGER  :: step_number = 0         ! The point-kinetic step number.
   REAL(dp) :: time_at_eos = 0.0_dp    ! The time at the end of step. (s)
   INTEGER  :: time_reject_count   = 0 ! Number of times the step size was rejected during step calculation.
   INTEGER  :: implicit_iter_count = 0 ! Implicit iterations performed.
   REAL(dp) :: omega = 0               ! Frequency from freq. transformation
   REAL(dp) :: omega_rt = 0.0_dp       ! Omega relative diff at the end of step (used in the conv. check)
   REAL(dp) :: n_eos_rt = 0.0_dp       ! Relative diff of the end of step neutron density (used in the conv. check)
END TYPE step_size_store
!
!  Reactor kinetic variables
!
   INTEGER,  SAVE, PROTECTED :: pc_group = 6 ! Number of delayed neutron precursor families (groups). (no unit)
   REAL(dp), SAVE, PROTECTED :: gen_time     ! Mean neutron generation time. (s)
   REAL(dp), SAVE, PROTECTED :: beta         ! Total delayed neutron fraction (sum over beta_l). (no unit)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: lambda_l(:)  ! Delayed neutron precusor decay constants. (1/s)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: beta_l(:)    ! Delayed neutron fraction. (no unit)
!
!  Solution control variables
!
   INTEGER,  SAVE, PROTECTED :: n_solution     = N_NO_FREQUENCY_TRANS ! Neutron density equation solution option.
   INTEGER,  SAVE, PROTECTED :: pre_approx     = C_FULLY_IMPLICIT ! Precursor concentration equation solution option.
   INTEGER,  SAVE, PROTECTED :: max_n_eos_iter = 25 ! Default max number of non-linear iterations on neutron dens.
   REAL(dp), SAVE, PROTECTED :: epsilon_n_eos  = 1.0e-7_dp ! Absolute neutron density convergence criteria.
   INTEGER,  SAVE, PROTECTED :: omega_option   = OMEGA_YES ! Determines how frequency is updated.
   INTEGER,  SAVE, PROTECTED :: converg_opt    = CONV_N    ! Conferge on the neutron density by default. 
!
!  Time control variables
!
   REAL(dp), SAVE, PROTECTED :: end_time       = 0.0_dp     ! Total desired simulation time. (s)
   REAL(dp), SAVE, PROTECTED :: tol_rel        = 1.0e-6_dp  ! Maximum accepted relative local truncation error. 
   REAL(dp), SAVE, PROTECTED :: min_ratio      = 0.20_dp    ! Minimum allowed ratio of the old/new time-step size.
   REAL(dp), SAVE, PROTECTED :: max_ratio      = 10.0_dp    ! Maximum allowed ratio of the old/new time-step size.
   REAL(dp), SAVE, PROTECTED :: edit_delta_t   = 0.1_dp     ! Time between printouts to file.(s)
   REAL(dp), SAVE, PROTECTED :: min_delta_t    = 1.0e-08_dp ! The minimum allowed time-step. (s)
   INTEGER,  SAVE, PROTECTED :: dt_option      = DT_ADAPT   ! Time control option. Adaptive time step is defualt
   REAL(dp), SAVE, PROTECTED :: dt_user        = 0.0_dp     ! User defiend constant time step. to be used with DT_CONST.
   INTEGER,  SAVE, PROTECTED :: time_reset_max = 3          ! The maximum number of time steps rejections pet step.
   REAL(dp), SAVE, PROTECTED :: exec_time      = 0.0_dp     ! Main loop execution time.
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: edit_table(:)  ! Table to store user-requested time edits.
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: knot_table(:)  ! Table to store all the time knots(combines edit_table,&
!                                                             regular edit_delta_t intervals, adn time_table_rho)
!
!  Reactivity variables
!
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: a_table(:)    ! Coefficients a for reactivity eq. Rho(t)=a(t)*t+b(t)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: b_table(:)    ! Coefficients b for reactivity eq. Rho(t)=a(t)*t+b(t)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: time_table_rho(:)   ! Table of time entries for a_table and b_table.
   CHARACTER(LEN=1), SAVE, PROTECTED      :: type_of_input = 'S' ! Type of reactivity input
!
!  Data storage linked list head and tail and some support variables
!
   TYPE(solution_storage), POINTER, SAVE, PROTECTED :: sol_list_head      ! Pointer to the first entry in the list.
   TYPE(solution_storage), POINTER, SAVE, PROTECTED :: sol_list_tail      ! Pointer to the last  entry in the list.
   INTEGER,                         SAVE, PROTECTED :: sol_list_count = 0 ! Number of entries in the linked list.
!
   TYPE(step_size_store), POINTER, SAVE, PROTECTED :: ss_list_head ! Pointer to the first entry in the list.
   TYPE(step_size_store), POINTER, SAVE, PROTECTED :: ss_list_top  ! Pointer the latest entry that was not yet fetched.
   TYPE(Step_size_store), POINTER, SAVE, PROTECTED :: ss_list_tail ! Pointer to the last  entry in the list.
   INTEGER,                        SAVE, PROTECTED :: ss_list_count = 0 ! Number of entries in the linked list.
!
!  Arrays for post-processing of the resuts (fs - final solution)
!
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: fs_time_stamps(:)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: fs_n(:)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: fs_rho(:)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: fs_c(:,:)
!
CONTAINS
!
!=======================================================================================================================
!
!  ************************
!  * set_data_real_single *
!  ************************
!
!  Purpose: Subroutine to change the value of a protected real variable (not an array!) in this module.
!  The following real variables can be modified:
!  gen_time, beta, lambda_l, beta_l, end_time, tol, edit_delta_t, min_delta_t, a_table, b_table, time_table_rho
!
!=======================================================================================================================
SUBROUTINE set_data_real_single(&
   name,                        & ! IN
   new_value,                   & ! IN
   idx)                           ! IN
!
!  Arguments
!
   CHARACTER(*)      , INTENT(IN) :: name      ! Variable name
   REAL(dp)          , INTENT(IN) :: new_value ! New value
   INTEGER , OPTIONAL, INTENT(IN) :: idx       ! Index of an entry in an array to be changed.
!
!
   SELECT CASE (name)
   CASE("gen_time")
      gen_time = new_value
   CASE("beta")
      beta = new_value
   CASE("lambda_l")
      IF ( PRESENT(idx) ) THEN
         lambda_l(idx) = new_value
      ELSE
         WRITE(*,*) 'An array index is needed to set a new value for lambda_l.'
      ENDIF
   CASE("beta_l")
      IF ( PRESENT(idx) ) THEN
         beta_l(idx) = new_value
      ELSE
         WRITE(*,*) 'An array index is needed to set a new value for beta_l.'
      ENDIF
   CASE("end_time")
      end_time = new_value
   CASE("exec_time")
      exec_time = new_value
   CASE("tol_rel")
      tol_rel = new_value
   CASE("edit_delta_t")
      edit_delta_t = new_value
   CASE("min_delta_t")
      min_delta_t = new_value
   CASE("epsilon_n_eos")
      epsilon_n_eos = new_value
   CASE("a_table")
      IF ( PRESENT(idx) ) THEN
         a_table(idx) = new_value
      ELSE 
         WRITE(*,*) 'An array index is needed to set a new value for a_table.'
      ENDIF
   CASE("b_table")
      IF ( PRESENT(idx) ) THEN
         b_table(idx) = new_value
      ELSE
         WRITE(*,*) 'An array index is needed to set a new value for b_table.'
      ENDIF
   CASE("time_table_rho")
      IF ( PRESENT(idx) ) THEN
         time_table_rho(idx) = new_value
      ELSE
         WRITE(*,*) 'An array index is needed to set a new value for time_table_rho.'
      ENDIF
   CASE("edit_table")
      IF ( PRESENT(idx) ) THEN
         edit_table(idx) = new_value
      ELSE
         WRITE(*,*) 'An array index is needed to set a new value for edit_table.'
      ENDIF
   END SELECT
!
END SUBROUTINE set_data_real_single
!=======================================================================================================================
!
!  ***********************
!  * set_data_real_array *
!  ***********************
!
!  Purpose: Subroutine to change the value of a protected real array in this module.
!  The following real variables can be modified:
!  lambda_l, beta_l, a_table, b_table, time_table
!
!=======================================================================================================================
SUBROUTINE set_data_real_array(&
   name,                       & ! IN
   new_value)                    ! IN
!
!  Arguments
!
   CHARACTER(*) , INTENT(IN) :: name         ! Variable name
   REAL(dp)     , INTENT(IN) :: new_value(:) ! New value
!
!
!
   SELECT CASE (name)
   CASE("lambda_l")
      lambda_l = new_value
   CASE("beta_l")
      beta_l = new_value
   CASE("a_table")
      a_table = new_value
   CASE("b_table")
      b_table = new_value
   CASE("time_table_rho")
      time_table_rho = new_value
   END SELECT
!
END SUBROUTINE set_data_real_array
!=======================================================================================================================
!
!  ****************
!  * set_data_int *
!  ****************
!
!  Purpose: Subroutine to change the value of a protected integer variable in this module.
!  The following integer variables can be modified:
!  pc_group, n_solution, max_n_eos_iter, pre_approx, time_reset_max
!
!=======================================================================================================================
SUBROUTINE set_data_int(&
   name,                & ! IN
   new_value)             ! IN
!
!  Arguments
!
   CHARACTER(*) , INTENT(IN) :: name      ! Variable name
   INTEGER      , INTENT(IN) :: new_value ! New value
!
!
!
   SELECT CASE (name)
   CASE("pc_group")
      pc_group = new_value
   CASE("n_solution")
      n_solution = new_value
   CASE("max_n_eos_iter")
      max_n_eos_iter = new_value
   CASE("pre_approx")
      pre_approx = new_value
   CASE("time_reset_max")
      time_reset_max = new_value
   CASE("omega_option")
      omega_option = new_value
   CASE("converg_opt")
      converg_opt = new_value
   END SELECT
!
END SUBROUTINE set_data_int
!=======================================================================================================================
!
!  ****************
!  * set_data_chr *
!  ****************
!
!  Purpose: Subroutine to change the value of a protected character variable in this module.
!  The following character variables can be modified:
!  type_of_input
!
!=======================================================================================================================
SUBROUTINE set_data_chr(&
   name,                & ! IN
   new_value)             ! IN
!
!  Arguments
!
   CHARACTER(*), INTENT(IN) :: name      ! Variable name
   CHARACTER(*), INTENT(IN) :: new_value ! New value
!
!
!
   SELECT CASE (name)
   CASE("type_of_input")
      IF ( LEN(new_value) > 1) THEN
         WRITE(*,*) 'type_of_input is expected to be a single character'
      ELSE
         type_of_input = new_value
      ENDIF
   END SELECT
!
END SUBROUTINE set_data_chr
!=======================================================================================================================
!
!  ******************
!  * set_dt_options *
!  ******************
!
!  Purpose: Subroutine that sets up the time step control options based on the user input.
!           Either adaptive time step control is used (no changes from the defualt values are required)
!           or the constant time step with value passed on by the user.
!
!=======================================================================================================================
SUBROUTINE set_dt_options(&
   delta_t_in)              ! IN
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: delta_t_in ! Desired constnat time step size (from user input)
!
   IF (delta_t_in > 0.0_dp) THEN
      dt_option = DT_CONST
      dt_user = delta_t_in
   ENDIF
!
END SUBROUTINE set_dt_options
!=======================================================================================================================
!
!  ***************
!  * allocate_1d *
!  ***************
!
!  Purpose: Subroutine to allocate 1D protected arrays in this module.
!  The following arrays can be allocated:
!  lambda_l, beta_l, a_table, b_table, time_table
!
!=======================================================================================================================
SUBROUTINE allocate_1d(&
   name,               & ! IN
   size)                 ! IN
!
!  Arguments
!
   CHARACTER(*), INTENT(IN) :: name ! Array name
   INTEGER     , INTENT(IN) :: size ! Size of the allocatable array
!
!
!
   SELECT CASE (name)
   CASE("lambda_l")
      ALLOCATE(lambda_l(size))
   CASE("beta_l")
      ALLOCATE(beta_l(size))  
   CASE("a_table")
      ALLOCATE(a_table(size))
   CASE("b_table")
      ALLOCATE(b_table(size))
   CASE("time_table_rho")
      ALLOCATE(time_table_rho(size))
   CASE("edit_table")
      ALLOCATE(edit_table(size))
   END SELECT
!
END SUBROUTINE allocate_1d
!=======================================================================================================================
!
!  *****************
!  * deallocate_1d *
!  *****************
!
!  Purpose: Subroutine to deallocate_1d 1D protected arrays in this module.
!  The following arrays can be deallocated:
!  lambda_l, beta_l, a_table, b_table, time_table_rho
!
!=======================================================================================================================
SUBROUTINE deallocate_1d(&
   name)                   ! IN
!
!  Arguments
!
   CHARACTER(*), INTENT(IN) :: name ! Array name
!
!
!
   SELECT CASE (name)
   CASE("lambda_l")
      DEALLOCATE(lambda_l)
   CASE("beta_l")
      DEALLOCATE(beta_l)
   CASE("a_table")
      DEALLOCATE(a_table)
   CASE("b_table")
      DEALLOCATE(b_table)
   CASE("time_table_rho")
      DEALLOCATE(time_table_rho)
   CASE("edit_table")
      DEALLOCATE(edit_table)
   END SELECT
!
END SUBROUTINE deallocate_1d
!=======================================================================================================================
!
!  ************************
!  * solution_storage_ini *
!  ************************
!
!  Purpose: Initialize the solution linked list.
!
!=======================================================================================================================
SUBROUTINE solution_storage_ini(&
   IN_solution,                 & ! IN
   IN_time_reject_count,        & ! IN
   IN_implicit_iter_count,      & ! IN
   IN_step_number,              & ! IN
   IN_time_at_eos)                ! IN
!
!  Arguments
!
   TYPE(solution_at_step), INTENT(IN) :: IN_solution
   INTEGER,  INTENT(IN) :: IN_time_reject_count   ! Number of times the step size was rejected during step calcualtion.
   INTEGER,  INTENT(IN) :: IN_implicit_iter_count ! Number of times the implicit iteration loop was entered.
   INTEGER,  INTENT(IN) :: IN_step_number         ! A number of the point-kinetic step in the simulation.
   REAL(dp), INTENT(IN) :: IN_time_at_eos         ! The time at the end of step.
!
   IF (.NOT. ASSOCIATED(sol_list_head)) THEN
      ALLOCATE(sol_list_head) ! Allocate new pointer
      sol_list_tail => sol_list_head
      NULLIFY(sol_list_tail%p_next_entry)
!
      sol_list_tail%kinetics_data = IN_solution
      sol_list_tail%time_reject_count = IN_time_reject_count
      sol_list_tail%implicit_iter_count = IN_implicit_iter_count
      sol_list_tail%step_number = IN_step_number
      sol_list_tail%time_at_eos = IN_time_at_eos
!
   ELSE
      WRITE(*,*) "ERROR *** The first entry of the linked list is already allocated!"
   ENDIF
!
   sol_list_count = sol_list_count + 1
!
END SUBROUTINE solution_storage_ini
!=======================================================================================================================
!
!  ************************
!  * solution_storage_add *
!  ************************
!
!  Purpose: Add a solution to the linked list of solutions.
!
!=======================================================================================================================
SUBROUTINE solution_storage_add(&
   IN_solution,                 & ! IN
   IN_time_reject_count,        & ! IN
   IN_implicit_iter_count,      & ! IN
   IN_step_number,              & ! IN
   IN_time_at_eos)                ! IN
!
!  Arguments
!
   TYPE(solution_at_step), INTENT(IN) :: IN_solution
   INTEGER,  INTENT(IN) :: IN_time_reject_count   ! Number of times the step size was rejected during step calcualtion.
   INTEGER,  INTENT(IN) :: IN_implicit_iter_count ! Number of times the implicit iteration loop was entered.
   INTEGER,  INTENT(IN) :: IN_step_number         ! A number of the point-kinetic step in the simulation.
   REAL(dp), INTENT(IN) :: IN_time_at_eos         ! The time at the end of step.
!
!  Calculations
!
   ALLOCATE(sol_list_tail%p_next_entry)        ! Allocate new pointer
   sol_list_tail => sol_list_tail%p_next_entry ! Last points to the new pointer
   NULLIFY(sol_list_tail%p_next_entry)         ! Nullify the pointer in the now last entry in the list
!
   sol_list_tail%kinetics_data = IN_solution
   sol_list_tail%time_reject_count = IN_time_reject_count
   sol_list_tail%implicit_iter_count = IN_implicit_iter_count
   sol_list_tail%step_number = IN_step_number
   sol_list_tail%time_at_eos = IN_time_at_eos
!
   sol_list_count = sol_list_count + 1
!
END SUBROUTINE solution_storage_add
!=======================================================================================================================
!
!  ****************************
!  * solution_storage_process *
!  ****************************
!
!  Purpose: Moves data from the solution storage into seprate 1D arrays.
!
!=======================================================================================================================
SUBROUTINE solution_storage_process()
!
!  Locals
!
   TYPE(solution_storage), POINTER :: p_entry ! Points to an (e)ntry i the linked list.
   INTEGER :: idx, idy ! Indecies in the final solution arrays.
!
!  Allocate memory for arrays holding the final solution (similariites to histrical events are coincidental).
!  -----------------------------------------------
   ALLOCATE(fs_time_stamps(sol_list_count))
   ALLOCATE(fs_n(sol_list_count))
   ALLOCATE(fs_rho(sol_list_count))
   ALLOCATE(fs_c(sol_list_count,pc_group))
!
!  Copy the data from the linked list to the arrays
!  -----------------------------------------------
   p_entry => sol_list_head
   procces_linekd_list: DO idx = 1, sol_list_count
      IF (.NOT. ASSOCIATED(p_entry)) THEN
         EXIT ! Pointer no longer valid. The end of the list has been reached.
      ENDIF
!
      fs_time_stamps(idx) = p_entry%time_at_eos
      fs_n(idx) = p_entry%kinetics_data%n2p
      fs_rho(idx) = p_entry%kinetics_data%rho
!
      DO idy = 1, pc_group
         fs_c(idx,idy) = p_entry%kinetics_data%c(idy)
      END DO
!
      p_entry => p_entry%p_next_entry
!
   END DO procces_linekd_list
!
END SUBROUTINE solution_storage_process
!=======================================================================================================================
!
!  ************************
!  * step_size_store_ini *
!  ************************
!
!  Purpose: Initialize the time step size linked list.
!
!=======================================================================================================================
SUBROUTINE step_size_store_ini()
!
   IF (.NOT. ASSOCIATED(ss_list_head)) THEN
      ALLOCATE(ss_list_head) ! Allocate new pointer
      ss_list_tail => ss_list_head
      ss_list_top => ss_list_head
      NULLIFY(ss_list_tail%p_next_entry)
!
   ELSE
      WRITE(*,*) "ERROR *** The first entry of the linked list is already allocated!"
   ENDIF
!
   ss_list_count = ss_list_count + 1
!
END SUBROUTINE step_size_store_ini
!=======================================================================================================================
!
!  ***********************
!  * step_size_store_add *
!  ***********************
!
!  Purpose: Add a solution to the linked list of solutions.
!
!=======================================================================================================================
SUBROUTINE step_size_store_add(&
   IN_step_number,             & ! IN
   IN_time_reject_count,       & ! IN
   IN_implicit_iter_count,     & ! IN
   IN_time_at_eos,             & ! IN
   IN_omega,                   & ! IN
   IN_omega_rt,                & ! IN 
   IN_n_eos_rt)                  ! IN
!
!  Arguments
!
   INTEGER , INTENT(IN) :: IN_step_number ! A number of the point-kinetic step in the simulation.
   REAL(dp), INTENT(IN) :: IN_time_at_eos ! The time at the end of step.
   INTEGER , INTENT(IN) :: IN_time_reject_count   ! Number of times the time step was rejected due to local error.
   INTEGER , INTENT(IN) :: IN_implicit_iter_count ! Number of iterations it took to reach convergence for a given step.
   REAL(dp), INTENT(IN) :: IN_omega               ! Frequency
   REAL(dp), INTENT(IN) :: IN_omega_rt ! Omega relative diff at the end of step (used in the conv. check)
   REAL(dp), INTENT(IN) :: IN_n_eos_rt ! Relative diff of the end of step neutron density (used in the conv. check)
!
!  Calculations
!
   ALLOCATE(ss_list_tail%p_next_entry)        ! Allocate new pointer
   ss_list_tail => ss_list_tail%p_next_entry ! Last points to the new pointer
   NULLIFY(ss_list_tail%p_next_entry)         ! Nullify the pointer in the now last entry in the list
!
   ss_list_tail%step_number = IN_step_number
   ss_list_tail%time_at_eos = IN_time_at_eos
   ss_list_tail%time_reject_count = IN_time_reject_count
   ss_list_tail%implicit_iter_count = IN_implicit_iter_count
   ss_list_tail%omega = IN_omega
   ss_list_tail%omega_rt = IN_omega_rt
   ss_list_tail%n_eos_rt = IN_n_eos_rt
!
   ss_list_count = ss_list_count + 1
!
END SUBROUTINE step_size_store_add
!=======================================================================================================================
!
!  *******************
!  * step_size_fetch *
!  *******************
!
!  Purpose: Recovers an entry from step size linked list.
!
!=======================================================================================================================
SUBROUTINE step_size_fetch(&
   step_number,            & ! OUT
   time_at_eos,            & ! OUT
   time_reject_count,      & ! OUT
   implicit_iter_count,    & ! OUT
   next_available,         & ! OUT
   omega,                  & ! OUT
   omega_rt,               & ! OUT
   n_eos_rt)                 ! OUT
!
!  Arguments
!
   INTEGER , INTENT(OUT) :: step_number
   REAL(dp), INTENT(OUT) :: time_at_eos
   INTEGER , INTENT(OUT) :: time_reject_count
   INTEGER , INTENT(OUT) :: implicit_iter_count
   LOGICAL , INTENT(OUT) :: next_available
   REAL(dp), INTENT(OUT) :: omega
   REAL(dp), INTENT(OUT) :: omega_rt
   REAL(dp), INTENT(OUT) :: n_eos_rt
! 
!  Copy the data from the linked list to the arrays
!  -----------------------------------------------
   step_number = ss_list_top%step_number
   time_at_eos = ss_list_top%time_at_eos
   time_reject_count   = ss_list_top%time_reject_count 
   implicit_iter_count = ss_list_top%implicit_iter_count
   omega    = ss_list_top%omega
   omega_rt = ss_list_top%omega_rt
   n_eos_rt = ss_list_top%n_eos_rt
!
   IF (.NOT. ASSOCIATED(ss_list_top%p_next_entry)) THEN
      next_available = .FALSE.
   ELSE
      next_available = .TRUE.
   ENDIF
!
   ss_list_top => ss_list_top%p_next_entry
!
END SUBROUTINE step_size_fetch
!=======================================================================================================================
!
!  *******************
!  * make_knot_table *
!  *******************
!
!  Purpose: Creates a knot table.
!
!=======================================================================================================================
SUBROUTINE make_knot_table()
!
!  Locals
!
   INTEGER :: size_temp = 0
   INTEGER :: size_edit
   INTEGER :: size_rho
   INTEGER :: size_all
   INTEGER :: idx ! Index variable     
   REAL(dp), ALLOCATABLE :: temp_a(:) ! Temporary array 1
 !  REAL(dp), ALLOCATABLE :: temp_a2(:) ! Temporary array 2
!
!  Calculate nubmer of entries needed between 
!   0 secodns end_time with edit_delta_t steps.
!
!  If the edit_delta_t was set to 0, it means no regular edit steps were requested.
   IF (edit_delta_t > 0.0_dp) THEN
      size_temp = FLOOR(end_time/edit_delta_t)+1
      ! If modulo is greater than 0, increase number of steps by 1 to account for the end_time
      IF (MOD(end_time,edit_delta_t) > 0.0_dp) THEN
         size_temp = size_temp+1
      ENDIF
   ELSE
      size_temp = 0
   ENDIF
!
!  number_of_steps now has all the 0 and end_time points accounted for
!  If not edit_time_t was defined simply add 0 and end_time to the temp_a
!
   IF (size_temp > 0) THEN
      ALLOCATE(temp_a(size_temp))
      DO idx = 1, size_temp-1
         temp_a(idx) = (idx-1)*edit_delta_t
      END DO
   ELSE
      size_temp = 2
      ALLOCATE(temp_a(size_temp))
   ENDIF
!
   temp_a(1) = 0.0_dp
   temp_a(size_temp) = end_time
!
!  Get sizes of all component arrays
!  -----------------------------------------------
   size_temp = SIZE(temp_a)
   size_rho = SIZE(time_table_rho)
!
   IF (ALLOCATED(edit_table)) THEN
      size_edit = SIZE(edit_table)
      size_all = size_temp + size_rho + size_edit
      ALLOCATE(knot_table(size_all))
      ! Populate the temporary table 2
      knot_table(1:size_temp) = temp_a
      knot_table((size_temp+1):(size_temp+size_rho)) = time_table_rho
      knot_table((size_temp+size_rho+1):size_all) = edit_table
!
   ELSE
      WRITE(*,*) "No time knots input file was defined."
      size_all = size_temp + size_rho
      ALLOCATE(knot_table(size_all))
      ! Populate the temporary table 2
      knot_table(1:size_temp) = temp_a
      knot_table((size_temp+1):size_all) = time_table_rho
   ENDIF
!
!  Mark the duplicate entires in the table.
!  Sort the temporary table and count duplicate entries (entries within one minimum step size from each other)
!
   CALL pk_hpsort_eps(knot_table)
   CALL pk_remove_duplicates(knot_table)
!
END SUBROUTINE make_knot_table
!=======================================================================================================================
END MODULE pk_data_types
!=======================================================================================================================