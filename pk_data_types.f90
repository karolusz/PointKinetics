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
!  CALL solution_storage_add()  to add an entry to the linked list.
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
   USE pk_enumeration, ONLY: C_FULLY_IMPLICIT, C_CONSTANT_APPROX, C_LINEAR_APPROX , C_EXPONENT_APPROX, &
                             N_FULLY_IMPLICIT, N_FREQUENCY_TRANS
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
   PUBLIC :: pc_group      , gen_time      , beta         , lambda_l     , beta_l
   PUBLIC :: n_solution    , pre_approx    , epsilon_n_eos, max_n_eos_iter
   PUBLIC :: end_time      , time_reset_max, eps_rel      , edit_delta_t , min_ratio , max_ratio, min_delta_t
   PUBLIC :: a_table       , b_table       , time_table   , type_of_input
   PUBLIC :: sol_list_head , sol_list_tail
   PUBLIC :: fs_time_stamps, fs_n          , fs_rho    , fs_c
!
!  Public subroutines
!  ------------------
   PUBLIC :: set_data
   PUBLIC :: allocate_1d, deallocate_1d
   PUBLIC :: solution_storage_ini, solution_storage_add, solution_storage_process
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
   REAL(dp) :: n2p = 1.0_dp ! Neutron (relative) density. Second-order method. (no unit)
   REAL(dp) :: n1p = 1.0_dp ! Neutron (relative) density. First-order method. (no unit)
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
   INTEGER,  SAVE, PROTECTED :: n_solution     = N_FULLY_IMPLICIT ! Neutron density equation solution option.
   INTEGER,  SAVE, PROTECTED :: pre_approx     = C_FULLY_IMPLICIT ! Precursor concentration equation solution option.
   INTEGER,  SAVE, PROTECTED :: max_n_eos_iter = 25 ! Default max number of non-linear iterations on neutron dens.
   REAL(dp), SAVE, PROTECTED :: epsilon_n_eos  = 1.0e-10_dp ! Absolute neutron density convergence criteria.
!
!  Time control variables
!
   REAL(dp), SAVE, PROTECTED :: end_time       = 0.0_dp     ! Total desired simulation time. (s)
   REAL(dp), SAVE, PROTECTED :: eps_rel        = 1.0e-08_dp ! Maximum accepted relative local truncation error. 
   REAL(dp), SAVE, PROTECTED :: min_ratio      = 0.5_dp     ! Minimum allowed ratio of the old/new time-step size.
   REAL(dp), SAVE, PROTECTED :: max_ratio      = 2.0_dp     ! Maximum allowed ratio of the old/new time-step size.
   REAL(dp), SAVE, PROTECTED :: edit_delta_t   = 0.1_dp     ! Time between printouts to file.(s)
   REAL(dp), SAVE, PROTECTED :: min_delta_t    = 1.0e-08_dp ! The minimum allowed time-step. (s)
   INTEGER,  SAVE, PROTECTED :: time_reset_max = 3          ! The maximum number of time steps rejections pet step.
!
!  Reactivity variables
!
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: a_table(:)    ! Coefficients a for reactivity eq. Rho(t)=a(t)*t+b(t)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: b_table(:)    ! Coefficients b for reactivity eq. Rho(t)=a(t)*t+b(t)
   REAL(dp), SAVE, PROTECTED, ALLOCATABLE :: time_table(:) ! Table of time entries for a_table and b_table.
   CHARACTER(LEN=1), SAVE, PROTECTED      :: type_of_input = 'S' ! Type of reactivity input
!
!  Data storage linked list head and tail and some support variables
!
   TYPE(solution_storage), POINTER, SAVE, PROTECTED :: sol_list_head      ! Pointer to the first entry in the list.
   TYPE(solution_storage), POINTER, SAVE, PROTECTED :: sol_list_tail      ! Pointer to the last entry in the list.
   INTEGER,                         SAVE, PROTECTED :: sol_list_count = 0 ! Number of entries in the linked list.
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
!  gen_time, beta, lambda_l, beta_l, end_time, eps_rel, edit_delta_t, min_delta_t, a_table, b_table, time_table
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
   CASE("eps_rel")
      eps_rel = new_value
   CASE("edit_delta_t")
      min_ratio = new_value
   CASE("min_delta_t")
      max_ratio = new_value
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
   CASE("time_table")
      IF ( PRESENT(idx) ) THEN
         time_table(idx) = new_value
      ELSE
         WRITE(*,*) 'An array index is needed to set a new value for time_table.'
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
   CASE("time_table")
      time_table = new_value
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
   CASE("time_table")
      ALLOCATE(time_table(size))
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
!  lambda_l, beta_l, a_table, b_table, time_table
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
   CASE("time_table")
      DEALLOCATE(time_table)
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
END MODULE pk_data_types
!=======================================================================================================================