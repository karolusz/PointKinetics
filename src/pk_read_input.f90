!=======================================================================================================================
!
!  *****************
!  * pk_read_input *
!  *****************
!
!  Purpose:
!  --------
!  Module to process the input files to the point kinetics program.
!
!  Usage:
!  ------
!  CALL process_input()  to process the kinetic parameters input file and reactor features input file.
!  CALL newunit()        to inquire about the available unit number.
!  CALL rho_dollar()     to calculate the reactivity at a given time.
!
!=======================================================================================================================
MODULE pk_read_input
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       08/03/18   K. Luszczek              Original code. Aggregates functions previously in separate modules.
!
!=======================================================================================================================
   USE pk_kinds     , ONLY: dp
   USE pk_data_types, ONLY: set_data, allocate_1d
!
   IMPLICIT NONE 
!
   PRIVATE
!
!  Public subroutines
!  ------------------
   PUBLIC :: process_input ! Subroutine which process the input files
   PUBLIC :: newunit       ! Function returns the available unit
   PUBLIC :: rho_dollar    ! Function that calculates rho dollar at a given time step
!
!=======================================================================================================================
CONTAINS
!=======================================================================================================================
!
!  ***********
!  * newunit *
!  ***********
!
!  Purpose: Function returns an available unit number for file processing.
! 
!=======================================================================================================================
INTEGER FUNCTION newunit(unit)
!
!  Arguments
!
   INTEGER, INTENT(OUT), OPTIONAL :: unit
!
!  Locals
!
   INTEGER :: LUN_MIN = 10, LUN_MAX = 1000
   LOGICAL :: is_opened
   INTEGER :: lun
!
!
!  ---------------------------------
   newunit = -1
   DO lun = LUN_MIN, LUN_MAX
      INQUIRE(unit=lun,opened=is_opened)
      IF(.NOT. is_opened) THEN
         newunit = lun
         EXIT
      ENDIF
   END DO
   IF(PRESENT(unit)) unit=newunit
!
END FUNCTION newunit
!=======================================================================================================================
!
!  **************
!  * rho_dollar *
!  **************
!
!  Purpose: Function returns the reactivity at a given time step.
!
!=======================================================================================================================
REAL(dp) FUNCTION rho_dollar(&
   time_current,             & !  IN
   type_of_input,            & !  IN
   a_table,                  & !  IN
   b_table,                  & !  IN
   time_table,               & !  IN
   lambda)                     !  IN, OPTIONAL
!=======================================================================================================================
! Record of revisions:
!     Date       Programmer               Description of changes
!     -----------------------------------------------------------
!     13/03/19   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp)        , INTENT(IN)  :: time_current  ! Current time in the simulation [s]
   CHARACTER(LEN=1), INTENT(IN)  :: type_of_input ! Reactivity input type
   REAL(dp)        , INTENT(IN)  :: a_table(:)    ! Table of a coefficients for linear reactivity eq. Rho(t) = a*t + b
   REAL(dp)        , INTENT(IN)  :: b_table(:)    ! Table of b coefficients for linear reactivity eq. Rho(t) = a*t + b
   REAL(dp)        , INTENT(IN)  :: time_table(:) ! Table of time entries for a_table and b_table.
!
!  Optional arguments (needed for a sinusoidal reactivity case)
!
   REAL(dp), OPTIONAL, INTENT(IN) :: lambda       ! 1-group delayed neutron precursor decay const. (1/s)
!                                                   For sinusoidal reactivity (Barry's paper section B3.2)
!                                                   "A refined way of solving reactor point kinetics equations for
!                                                   imposed reactivity insertions."
!
!  Locals
!
   INTEGER :: idx    ! Loop index variable
!
!  Check which reactivity equation is to be used for the current time-step
!  a_table, b_table and time_table all have the same size
   DO idx = 1, SIZE(time_table)
      IF (time_current <= time_table(idx)) THEN
!        Select the type of reactivity input
         SELECT CASE (type_of_input)
!
         CASE ("I") ! Step reactivity
!        ----------------------------
            rho_dollar = b_table(idx) + a_table(idx) * time_current
!
         CASE ("R") ! Ramp reactivity
!        ----------------------------
            rho_dollar = b_table(idx) + a_table(idx) * time_current
!
         CASE ("Z") ! Zigzag insertion
!        ----------------------------
            rho_dollar = b_table(idx) + a_table(idx) * time_current
!         
         CASE ("S") ! Sinusoidal insertion
!        ----------------------------
            rho_dollar = a_table(idx) / (a_table(idx) + lambda * b_table(idx))       &
                         * sin( 4.0_dp * atan (1.0_dp) * time_current / b_table(idx))
!
         END SELECT
!
         EXIT  ! Exit loop. The loop should not run after the first successful IF
!
      ENDIF
   END DO
END FUNCTION rho_dollar
!=======================================================================================================================
!
!  *****************
!  * process_input *
!  *****************
!
!  Purpose: Process the input file.
! 
!=======================================================================================================================
SUBROUTINE process_input(&
   kin_param_file_path,  & ! IN
   reactivity_file_path, & ! IN
   ierror)                 ! OUT
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       11/03/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   CHARACTER(*) , INTENT(IN)  :: kin_param_file_path  ! Kinetics parameters input filepath
   CHARACTER(*) , INTENT(IN)  :: reactivity_file_path ! Reactivity input filepath
   INTEGER      , INTENT(OUT) :: ierror ! Error number
!
!  Process the reactivity input file
!  ---------------------------------
   CALL process_reactor_input(&
      filename = kin_param_file_path, & ! IN
      ierror   = ierror)                ! OUT
!
!  Process the reactor input file
!  ------------------------------
   CALL process_reactivity_input(&
      filename    = reactivity_file_path, & ! IN
      ierror      = ierror)                 ! OUT
!
END SUBROUTINE process_input
!=======================================================================================================================
!
!  *************************
!  * process_reactor_input *
!  *************************
!
!  Purpose: Subroutine processes the reactor input file.
!
!=======================================================================================================================
SUBROUTINE process_reactor_input(&
   filename,                     & ! IN
   ierror)                         ! OUT
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!       11/03/18   K. Luszczek              Modified to use derived data types
!
!=======================================================================================================================
! The structure of reactor.in file should be as follows:
! Line 1: integer informing about the number of precursor groups, neutron generation time
! Line 2: group 1 lambda_l, group 1 beta_l
! Line n: group n lambda_l, group n beta_l
! ...
!
!  Arguments
!
   CHARACTER(*) , INTENT(IN)  :: filename ! Kinetics parameters input filepath
   INTEGER      , INTENT(OUT) :: ierror   ! Error handling number
!
!  Locals
!
   INTEGER :: idx      ! Loop index variable
   INTEGER :: unit_pkp ! Point kinetic parameter input file unit number
   REAL(dp) :: gen_time   ! Mean neutron generation time: variable to read to from file
   REAL(dp) :: lambda_l_r ! Variable to read to from file
   REAL(dp) :: beta_l_r   ! Variable to read to from file 
   REAL(dp) :: beta_r     ! Variable to read to from file
   INTEGER :: pc_group_r  ! Number of pre-cursor groups
!
!  Open the reactivity file
!  ------------------------
   unit_pkp = newunit()
!
   OPEN (UNIT = unit_pkp, FILE= filename, STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierror)
   IF (ierror /= 0) THEN
      WRITE(*,*) 'Error occurred when opening point-kinetics parameter file. IOSTAT error: ', ierror
      RETURN
   ENDIF   
!
!  Read number of lines and time step
!  ----------------------------------
   READ(unit_pkp,*, IOSTAT = ierror) pc_group_r, gen_time
   IF (ierror /= 0) THEN
      WRITE(*,*) 'Error occurred when reading point-kinetics parameter file. IOSTAT error: ', ierror
      RETURN
   ENDIF 
!
   CALL set_data('pc_group' ,pc_group_r)
   CALL set_data('gen_time' ,gen_time)
   CALL allocate_1d('lambda_l' ,pc_group_r)
   CALL allocate_1d('beta_l' ,pc_group_r)
!
   beta_r = 0.0_dp
!
   DO idx = 1, pc_group_r
      READ(unit_pkp,*, IOSTAT = ierror) lambda_l_r, beta_l_r
      IF (ierror /= 0) THEN
         WRITE(*,*) 'Error occurred when reading point-kinetics parameter file. IOSTAT error: ', ierror
         RETURN
!
      ELSE
!
         beta_r = beta_r + beta_l_r
         CALL set_data('lambda_l', lambda_l_r, idx)
         CALL set_data('beta_l', beta_l_r , idx)
!
      ENDIF
   END DO
! 
   CLOSE (UNIT = unit_pkp)
!
   CALL set_data('beta',beta_r)
!
END SUBROUTINE process_reactor_input
!=======================================================================================================================
!
!  ****************************
!  * process_reactivity_input *
!  ****************************
!
!  Purpose: Process the reactivity input file.
!
!=======================================================================================================================
SUBROUTINE process_reactivity_input(&
   filename,                        & ! IN
   ierror)                            ! OUT
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!       11/03/18   K. Luszczek              Modified to use derived data types
!
!=======================================================================================================================
!
!  Arguments
!
   CHARACTER(*) , INTENT(IN)  :: filename ! File path to the reactivity input
   INTEGER      , INTENT(OUT) :: ierror   ! Error handling number
!
!  Locals
!
   INTEGER  :: idx                   ! Loop index variable
   REAL(dp) :: time                  ! Total simulation time
   INTEGER  :: lines                 ! Number of lines to be read from the input file
   REAL(dp) :: edit_delta_t          ! Time step
   INTEGER  :: unit_ri               ! Reactivity input file unit number
   REAL(dp) :: a_table_r             ! Variable to store a_table entry from the input file.
   REAL(dp) :: b_table_r             ! Variable to store b_table entry from the input file.
   REAL(dp) :: time_table_r          ! Variable to store time_table entry from the input file.
   INTEGER  :: n_solution_in         ! Neutron density equation solution option
   INTEGER  :: pre_approx_in         ! Pre-cursor concentration equation solution option
   REAL(dp) :: epsilon_n_eos_in      ! Neutron density convergence criteria
   CHARACTER(LEN=1) :: type_of_input ! Type of reactivity input
!
! Open the reactivity file
!  ------------------------
   unit_ri = newunit()
!
   OPEN (UNIT = unit_ri, FILE = filename, STATUS = 'OLD', ACTION = 'READ', IOSTAT = ierror) ! Open the reactivity file
   IF (ierror /= 0) THEN
      WRITE(*,*) 'Error occurred when opening reactivity file. IOSTAT error: ', ierror
      RETURN
   ENDIF     
!
   READ(unit_ri, *, IOSTAT = ierror) type_of_input, lines, edit_delta_t, n_solution_in, pre_approx_in, epsilon_n_eos_in
   IF (ierror /= 0) THEN
      WRITE(*,*) 'Error occurred when reading reactivity file. IOSTAT error: ', ierror
      RETURN
   ENDIF 
!
   CALL set_data("edit_delta_t",edit_delta_t)
   CALL allocate_1d("a_table", lines)
   CALL allocate_1d("b_table", lines)
   CALL allocate_1d("time_table", lines)
!
   process_input_select: SELECT CASE (type_of_input)
!
   CASE ("I","R","Z") ! Step reactivity, ramp reactivity, zigzag insertion
!  ---------------------------- 
      DO idx = 1, lines
         READ(unit_ri,*, IOSTAT = ierror) b_table_r, a_table_r, time_table_r
         IF (ierror /= 0) THEN
            WRITE(*,*) 'Error occurred when reading reactivity file. IOSTAT error: ', ierror
            RETURN
         ELSE
            CALL set_data("b_table", b_table_r, idx)
            CALL set_data("a_table", a_table_r, idx)
            CALL set_data("time_table", time_table_r, idx)
         ENDIF         
      END DO   
      time = time_table_r !  The last entry from the file defines the total simulation time
!
   CASE ("S") ! Sinusoidal insertion
!  ---------------------------- 
      READ(unit_ri,*, IOSTAT = ierror) a_table_r, b_table_r, time_table_r
      IF (ierror /= 0) THEN
         WRITE(*,*) 'Error occurred when reading reactivity file. IOSTAT error: ', ierror
         RETURN
      ELSE
         CALL set_data("a_table", a_table_r, 1)
         CALL set_data("b_table", b_table_r, 1)
         CALL set_data("time_table", time_table_r, 1)
      ENDIF  
      time = time_table_r !  The last entry from the file defines the total simulation time
!
   CASE DEFAULT ! no case recognized
!  ----------------------------
      WRITE(*,*) ('Error: reactivity case not recognized.')
      ierror = 1
!
   END SELECT process_input_select
   CLOSE (UNIT = unit_ri) 
!
!  Move loaded data to the data structure
!  ----------------------------------
   CALL set_data("end_time", time)
!
   CALL set_data("type_of_input", type_of_input)
!
   CALL set_data("n_solution", n_solution_in)
   CALL set_data("pre_approx", pre_approx_in)
   CALL set_data("epsilon_n_eos", epsilon_n_eos_in)
!
END SUBROUTINE process_reactivity_input
!=======================================================================================================================
END MODULE pk_read_input
!=======================================================================================================================
