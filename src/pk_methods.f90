!=======================================================================================================================
!
!  **************
!  * pk_methods *
!  **************
!
!  Purpose:
!  --------
!
!  This module calculates the neutron flux and delayed neutron precursor concentrations.
!
!  Usage:
!  ------
!
!  CALL point_kinetics_step   to perform one step in the simulation
!
!=======================================================================================================================
MODULE pk_methods
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!  Import subroutines
   USE pk_read_input  , ONLY: rho_dollar
   USE pk_write_output, ONLY: pk_write_to_file
!  Import variables
   USE pk_kinds       , ONLY: dp
   USE pk_enumeration , ONLY: C_FULLY_IMPLICIT    , C_CONSTANT_APPROX, C_LINEAR_APPROX , C_EXPONENT_APPROX, &
                              N_NO_FREQUENCY_TRANS, N_FREQUENCY_TRANS, OMEGA_NO        , OMEGA_YES        , &
                              CONV_N              , CONV_W
   USE pk_data_types  , ONLY: set_data            , solution_at_step , pc_group        , max_n_eos_iter   , &
                              n_solution          , pre_approx       , epsilon_n_eos   , omega_option     , &
                              converg_opt
!
   IMPLICIT NONE
!
   PRIVATE
!
!  Public subroutines
!  ------------------
   PUBLIC :: point_kinetics_step ! Subroutine that performs one point kinetic time step
   PUBLIC :: pk_methods_init     ! Initializes the this module.
!
!  Module variables
!  ------------------
   REAL(dp), SAVE :: omega       ! Neutron density frequency. Set to zero at first. Save between the time steps.
!
!=======================================================================================================================
CONTAINS
!=======================================================================================================================
!
!  ***********************
!  * point_kinetics_step *
!  ***********************
!
!  Purpose: calculate one time step. Returns updated n_eos and n_eos_p1 (calculated with 2nd and 1st order methods)
!           as well as c_eos.
! 
!=======================================================================================================================
SUBROUTINE point_kinetics_step(&
   gen_time,                   & ! IN (s)
   lambda_l,                   & ! IN (1/s)
   beta,                       & ! IN (-)
   beta_l,                     & ! IN (-)
   n_bos,                      & ! IN (1/cm^3)
   c_bos,                      & ! IN (1/cm^3)
   rho_eos,                    & ! IN (-)
   rho_bos,                    & ! IN (-)
   delta_t,                    & ! IN (s)
   n_eos,                      & ! OUT (1/cm^3)
   n_eos_p1,                   & ! OUT (1/cm^3)
   n_eos_tsc,                  & ! OUT (1/cm^3 or -)
   n_eos_p1_tsc,               & ! OUT (1/cm^3 or -)
   c_eos,                      & ! OUT (1/cm^3)
   inner_iter,                 & ! OUT (1/s)
   omega_out,                  & ! OUT (-)
   omega_rt,                   & ! OUT (-)
   n_eos_rt )                    ! OUT (-)
!=======================================================================================================================
! Record of revisions:
!       Date     Programmer  Description of change
!       ====     ==========  =====================
!       27/09/18 K. Luszczek Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN)  :: gen_time    ! Mean neutron generation time                     (s)
   REAL(dp), INTENT(IN)  :: lambda_l(:) ! Table to store precursor decay constants         (1/s)
   REAL(dp), INTENT(IN)  :: beta        ! Total delayed neutron fraction                   (no unit)
   REAL(dp), INTENT(IN)  :: beta_l(:)   ! Group-wise core delayed neutron fractions        (no unit)
   REAL(dp), INTENT(IN)  :: n_bos       ! Initial core average neutron density at t(i)     (1/cm^3)
   REAL(dp), INTENT(IN)  :: c_bos(:)    ! initial group-wise precursors' concentr. at t(i) (1/cm^3)
   REAL(dp), INTENT(IN)  :: rho_eos     ! Reactivity, at the end of step                   (-)
   REAL(dp), INTENT(IN)  :: rho_bos     ! Reactivity at the beginning of step              (-)
   REAL(dp), INTENT(IN)  :: delta_t     ! Time step size                                   (1/s)
!
   REAL(dp), INTENT(OUT) :: n_eos      ! Core average neutron density at EOS from the 2nd order method   (1/cm^3)
   REAL(dp), INTENT(OUT) :: n_eos_p1   ! Core average neutron density at EOS from the 1st order method   (1/cm^3)
   REAL(dp), INTENT(OUT) :: n_eos_tsc  ! Core average n. dens. at EOS (for time step control). 2nd order method.
   REAL(dp), INTENT(OUT) :: n_eos_p1_tsc ! Core average n. dens. at EOS (for time step control). 1st order method.
   REAL(dp), INTENT(OUT) :: c_eos(:)   ! Group-wise precursors' concentration at EOS                     (1/cm^3)
   REAL(dp), INTENT(OUT) :: omega_out  ! Frequecny                                                       (1/s)
   INTEGER,  INTENT(OUT) :: inner_iter ! Count of inner iterations.                                      (no unit)
   REAL(dp), INTENT(OUT) :: omega_rt   ! Relative difference between omega iterations (from conv. check))(-)
   REAL(dp), INTENT(OUT) :: n_eos_rt   ! Relative difference between iterations of n_eos                 (-)
!
!  Locals
! 
   REAL(dp) :: x_tilde(pc_group) ! c_eos equation coefficient
   REAL(dp) :: y_tilde(pc_group) ! c_eos equation coefficient
   REAL(dp) :: z_tilde(pc_group) ! c_eos equation coefficient 
!
   REAL(dp) :: theta ! Theta value detemines the order of the method. 1.0 -> 1st order, 0.5 -> 2nd order.
!
!  Calculations
!
! Set n_eos and c_eos equal to n_bos and c_bos
! --------------------------------------------
   c_eos(:) = c_bos(:)
   n_eos = n_bos
!
!  Select the fission source approximation
!  --------------------------------------------
   CALL INNER_pc_meth_select()
!
!  MAIN LOOP FOR N_EOS ITERATIONS
!  --------------------------------------------
   CALL INNER_n_eos_main_loop()
! 
!  Update c_eos (precursor concentration at the end-of-step)
!  --------------------------------------------
   CALL update_c_eos(&
      gen_time = gen_time, & ! IN 
      lambda_l = lambda_l, & ! IN
      beta_l   = beta_l  , & ! IN
      x_tilde  = x_tilde , & ! IN
      y_tilde  = y_tilde , & ! IN
      z_tilde  = z_tilde , & ! IN
      n_bos    = n_bos   , & ! IN
      n_eos    = n_eos   , & ! IN
      c_bos    = c_bos   , & ! IN
      c_eos    = c_eos)      ! OUT
! 
!  Calculate the n_eos using first order methods. 
!  This value will be used for local truncation error estimate.
!  --------------------------------------------
   CALL INNER_first_order_n()
!
!  Determine eos values for time step control
!  --------------------------------------------
   CALL INNER_time_step_control_variables()
!
!  --------------------------------------------
!  Internal subroutines
!  --------------------------------------------
!
   CONTAINS
!=======================================================================================================================
   SUBROUTINE INNER_pc_meth_select()
!=======================================================================================================================
!
!  The only two cases when the iterations are needed on EOS neutron density are when:
!  1. Exponential approxmation of the fission source is selected.
!  2. Frequency trasnformed solution method is selected.
!  Thus, implicit, constant, linear fission source approximations (those that do not depend on EOS neutron dens.)
!  are taken out of the loop to avoid unnecessary reculacation during iterations.
!
      pc_meth_select_outside_loop: SELECT CASE (pre_approx)
!
         CASE (C_FULLY_IMPLICIT) ! fully implicit Euler
!
            CALL pc_thimp_approx(&
               lambda_l = lambda_l, & ! IN (1/s)
               delta_t  = delta_t , & ! IN (s)
               x_tilde  = x_tilde , & ! OUT (-)
               y_tilde  = y_tilde , & ! OUT (-)
               z_tilde  = z_tilde)    ! OUT (-)
!
         CASE (C_CONSTANT_APPROX) ! constant approximation of the source in precursor eq.
!
            CALL pc_const_approx(&
               lambda_l = lambda_l, & ! IN (1/s)
               delta_t  = delta_t , & ! IN (s)
               x_tilde  = x_tilde , & ! OUT (-)
               y_tilde  = y_tilde , & ! OUT (-)
               z_tilde  = z_tilde)    ! OUT (-)
! 
         CASE (C_LINEAR_APPROX) ! linear approximation of the source in precursor eq.
!
            CALL pc_linea_approx(&
               lambda_l = lambda_l, & ! IN (1/s)
               delta_t  = delta_t , & ! IN (s)
               x_tilde  = x_tilde , & ! OUT (-)
               y_tilde  = y_tilde , & ! OUT (-)
               z_tilde  = z_tilde)    ! OUT (-)
!
         CASE (C_EXPONENT_APPROX) ! exponential approximation of the source in precursor eq.
!
            CONTINUE ! Will be handled inside the n_eos_loop if needed.
!
         CASE DEFAULT ! Issue a warning
!
            WRITE(*,*) ' ERROR ** No proper fission source approximation was seleced!.'
!
      END SELECT pc_meth_select_outside_loop
   END SUBROUTINE INNER_pc_meth_select
!=======================================================================================================================
   SUBROUTINE INNER_n_eos_main_loop()
!=======================================================================================================================
!
!     Locals
!
      INTEGER :: idx ! Loop index variable.
      REAL(dp) :: n_eos_i1           ! Previous iteration step neutron density         (1/cm^3)
      REAL(dp) :: c_eos_i1(pc_group) ! Previous iteration step precursor concentration (1/cm^3)
      REAL(dp) :: omega_i1           ! Previous iteration step omega (-)
!
!     MAIN LOOP FOR N_EOS ITERATIONS
!     -------------------------------
      theta = 0.5_dp
!
      loop_implicit: DO idx = 1, max_n_eos_iter
!
         n_eos_i1= n_eos        ! Save pre-iteration n_eos value.
         c_eos_i1(:) = c_eos(:) ! Save pre-iteration c_eos value.
         inner_iter = idx       ! Track number of iterations performed by the subroutine.
!
!     Check if the exponential fission source approximation was selected.
!     -------------------------------------------------------------------
      IF (pre_approx == C_EXPONENT_APPROX) THEN
!
         CALL pc_expon_approx(&
            lambda_l = lambda_l, & ! IN (1/s)
            delta_t  = delta_t , & ! IN (s)
            n_bos    = n_bos   , & ! IN (1/cm^3)
            n_eos    = n_eos   , & ! IN (1/cm^3)
            omega_p  = omega   , & ! IN (-)
            x_tilde  = x_tilde , & ! OUT (-)
            y_tilde  = y_tilde , & ! OUT (-)
            z_tilde  = z_tilde)    ! OUT (-)
!
      ENDIF
!
!     If n equation option is not frequency transformed: set omega to 0.
!                                             otheriwse: use actual omega
!     -------------------------------------------------------------------
      IF (n_solution == N_NO_FREQUENCY_TRANS) THEN
         CALL update_n_eos(&
            gen_time = gen_time, & ! IN (s)
            lambda_l = lambda_l, & ! IN (1/s)
            delta_t  = delta_t , & ! IN (s)
            rho_eos  = rho_eos , & ! IN (-)
            rho_bos  = rho_bos , & ! IN (-)
            beta     = beta    , & ! IN (-)
            beta_l   = beta_l  , & ! IN (-)
            x_tilde  = x_tilde , & ! IN (-)
            y_tilde  = y_tilde , & ! IN (-)
            z_tilde  = z_tilde , & ! IN (-)
            omega_n  = 0.0_dp  , & ! IN (-)
            theta    = theta   , & ! IN (-)
            c_bos    = c_bos   , & ! IN (-)
            n_bos    = n_bos   , & ! IN (1/cm^3)
            n_eos    = n_eos)      ! OUT (1/cm^3)
!
      ELSE
         CALL update_n_eos(&
            gen_time = gen_time, & ! IN (s)
            lambda_l = lambda_l, & ! IN (1/s)
            delta_t  = delta_t , & ! IN (s)
            rho_eos  = rho_eos , & ! IN (-)
            rho_bos  = rho_bos , & ! IN (-)
            beta     = beta    , & ! IN (-)
            beta_l   = beta_l  , & ! IN (-)
            x_tilde  = x_tilde , & ! IN (-)
            y_tilde  = y_tilde , & ! IN (-)
            z_tilde  = z_tilde , & ! IN (-)
            omega_n  = omega   , & ! IN (-)
            theta    = theta   , & ! IN (-)
            c_bos    = c_bos   , & ! IN (-)
            n_bos    = n_bos   , & ! IN (1/cm^3)
            n_eos    = n_eos)      ! OUT (1/cm^3)      
!
      ENDIF
! 
!     Check convergence of n_eos or if only one pass (fully implicit) was requested. if 'yes' to any, break the loop.
!
         110 FORMAT ( I2, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6 )
         WRITE(*,110) inner_iter, delta_t, n_bos, n_eos, omega, rho_bos, rho_eos
!
         IF (omega_option == OMEGA_YES) THEN
            omega_i1 = omega
            omega = log(n_eos/n_bos) / delta_t ! Update the omega value.
         ENDIF
!
         omega_rt = abs(omega - omega_i1)/omega
         n_eos_rt = abs(n_eos - n_eos_i1)/n_eos
!
         convergence_check: SELECT CASE (converg_opt)
!
            CASE (CONV_N)
!           
               IF ((n_eos_rt <= epsilon_n_eos) .OR. (idx == max_n_eos_iter)) THEN
                  omega_out = omega
                  EXIT loop_implicit
               ENDIF
!
            CASE (CONV_W)
!
               IF ((omega_rt <= epsilon_n_eos) .OR. (idx == max_n_eos_iter)) THEN
                  omega_out = omega
                  EXIT loop_implicit
               ENDIF
!
         END SELECT convergence_check
!
      END DO loop_implicit
   END SUBROUTINE INNER_n_eos_main_loop
!=======================================================================================================================
   SUBROUTINE INNER_first_order_n()
!=======================================================================================================================
!     Calculate the n_eos using the first order method (for local truncation error estimate later on).
!     Fully implicit method is assumed (no iterations needed). The omega is set back to 1.0 for methotds without 
!     frequency transformation.
!     ------------
      theta = 1.0_dp
      IF (n_solution == N_NO_FREQUENCY_TRANS) THEN
         omega = 1.0_dp
      ENDIF
!
      CALL update_n_eos(&
         gen_time = gen_time , & ! IN (s)
         lambda_l = lambda_l , & ! IN (1/s)
         delta_t  = delta_t  , & ! IN (s)
         rho_eos  = rho_eos  , & ! IN (-)
         rho_bos  = rho_bos  , & ! IN (-)
         beta     = beta     , & ! IN (-)
         beta_l   = beta_l   , & ! IN (-)
         x_tilde  = x_tilde  , & ! IN (-)
         y_tilde  = y_tilde  , & ! IN (-)
         z_tilde  = z_tilde  , & ! IN (-)
         omega_n  = omega    , & ! IN (-)
         theta    = theta    , & ! IN (-)
         c_bos    = c_bos    , & ! IN (-)
         n_bos    = n_bos    , & ! IN (1/cm^3)
         n_eos    = n_eos_p1)    ! OUT (1/cm^3)
!
   END SUBROUTINE INNER_first_order_n
!=======================================================================================================================
   SUBROUTINE INNER_time_step_control_variables()
!=======================================================================================================================
!     Sets the approptiate values for variables used in the time-step control. For non frequency transformed methods
!     these are equal to end of step neutron densiities. For frequncy transfromed methods these should be transformed
!     rather than abosulte values.
!     ------------
      IF (n_solution == N_NO_FREQUENCY_TRANS) THEN
         n_eos_tsc = n_eos
         n_eos_p1_tsc = n_eos_p1
      ELSE
         n_eos_tsc = n_eos * exp(-1.0_dp * omega * delta_t)
         n_eos_p1_tsc = n_eos_p1 * exp(-1.0_dp * omega * delta_t)
      ENDIF
!
   END SUBROUTINE INNER_time_step_control_variables

END SUBROUTINE point_kinetics_step
!=======================================================================================================================
!
!  *******************
!  * pc_thimp_approx *
!  *******************
!
!  Purpose: Subroutine to calculate x_tilde, y_tilde, and z_tilde for the fully implicit
!           backward Euler integration of the precursor equation.
!=======================================================================================================================
SUBROUTINE pc_thimp_approx(&
   lambda_l,               & ! IN (1/s)
   delta_t,                & ! IN (s)
   x_tilde,                & ! OUT (-)
   y_tilde,                & ! OUT (-)
   z_tilde)                  ! OUT (-)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Precursor decay constants (1/s)
   REAL(dp), INTENT(IN) :: delta_t     ! Time step size            (s)
!
   REAL(dp), INTENT(OUT) :: x_tilde(:) ! Factor in EOS delayed n. pre. density equation (-)
   REAL(dp), INTENT(OUT) :: y_tilde(:) ! Factor in EOS delayed n. pre. density equation (-)
   REAL(dp), INTENT(OUT) :: z_tilde(:) ! Factor in EOS delayed n. pre. density equation (-)
!
!  Locals
!
   INTEGER  :: idx                   ! Loop index variable
   REAL(dp) :: lambda_star(pc_group) ! Variable need for calculations
!
!  Calculations
!
   thimp_approx: DO idx = 1, pc_group
!
      lambda_star(idx) = lambda_l(idx) / (1.0_dp + delta_t * lambda_l(idx))
!
      x_tilde(idx) = lambda_star(idx) / lambda_l(idx)
      y_tilde(idx) = 0.0_dp
      z_tilde(idx) = delta_t * lambda_star(idx)
!
   END DO thimp_approx
END SUBROUTINE pc_thimp_approx
!=======================================================================================================================
!
!  *******************
!  * pc_const_approx *
!  *******************
!
!  Purpose: Subroutine to calculate x_tilde, y_tilde, and z_tilde the constant fission source approximation
!           in the analytical precursor group integration.
!=======================================================================================================================
SUBROUTINE pc_const_approx(&
   lambda_l,               & ! IN (1/s)
   delta_t,                & ! IN (s)
   x_tilde,                & ! OUT (-)
   y_tilde,                & ! OUT (-)
   z_tilde)                  ! OUT (-)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Precursor decay constants (1/s)
   REAL(dp), INTENT(IN) :: delta_t     ! Time step size            (1/s)
!
   REAL(dp), INTENT(OUT) :: x_tilde(:) ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(OUT) :: y_tilde(:) ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(OUT) :: z_tilde(:) ! Factor in EOS delayed n. pre. density equation
!
!  Locals
!
   INTEGER :: idx ! Loop index variable
!
!  Calculations
!
   const_approx: DO idx = 1, pc_group
!
      x_tilde(idx) = exp(-1.0_dp * lambda_l(idx) * delta_t)
      y_tilde(idx) = 1.0_dp - x_tilde(idx)
      z_tilde(idx) = 0.0_dp 
!
   END DO const_approx
END SUBROUTINE pc_const_approx
!=======================================================================================================================
!
!  *******************
!  * pc_linea_approx *
!  *******************
!
!  Purpose: Subroutine to calculatex_tilde, y_tilde, and z_tilde for the linear fission source approximation
!           in the analytical precursor group integration.
!=======================================================================================================================
SUBROUTINE pc_linea_approx(&
   lambda_l,               & ! IN (1/s)
   delta_t,                & ! IN (s)
   x_tilde,                & ! OUT (-)
   y_tilde,                & ! OUT (-)
   z_tilde)                  ! OUT (-)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Precursor decay constants (1/s)
   REAL(dp), INTENT(IN) :: delta_t     ! Time step size            (1/s)
!
   REAL(dp), INTENT(OUT) :: x_tilde(:) ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(OUT) :: y_tilde(:) ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(OUT) :: z_tilde(:) ! Factor in EOS delayed n. pre. density equation
!
!  Locals
!
   INTEGER :: idx ! Loop index variable
!
!  Calculations
!  
   linear_approx: DO idx = 1, pc_group
!
      x_tilde(idx) = exp(-1.0_dp * lambda_l(idx) * delta_t)
      z_tilde(idx) = 1.0_dp + (x_tilde(idx) - 1.0_dp) / (lambda_l(idx) * delta_t)
      y_tilde(idx) = 1.0_dp - x_tilde(idx) - z_tilde(idx)
!
   END DO linear_approx
END SUBROUTINE pc_linea_approx
!=======================================================================================================================
!
!  *******************
!  * pc_expon_approx *
!  *******************
!
!  Purpose: Subroutine to calculate x_tilde, y_tilde, and z_tilde for the exponential fission source approximation
!           in the analytical precursor group integration.
!=======================================================================================================================
SUBROUTINE pc_expon_approx(&
   lambda_l,               & ! IN (1/s)
   delta_t,                & ! IN (s)
   n_bos,                  & ! IN (1/cm^3)
   n_eos,                  & ! IN (1/cm^3)
   omega_p,                & ! IN
   x_tilde,                & ! OUT (-)
   y_tilde,                & ! OUT (-)
   z_tilde)                  ! OUT (-)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Table to store precursor decay constants  (1/s)
   REAL(dp), INTENT(IN) :: delta_t     ! Time step size (1/s)
   REAL(dp), INTENT(IN) :: n_bos       ! Beginning of step neutron dependence
   REAL(dp), INTENT(IN) :: n_eos       ! End of step neutron density
   REAL(dp), INTENT(IN) :: omega_p     ! Frequency
!
   REAL(dp), INTENT(OUT) :: x_tilde(:) ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(OUT) :: y_tilde(:) ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(OUT) :: z_tilde(:) ! Factor in EOS delayed n. pre. density equation
!
!  Locals
!
   INTEGER  :: idx               ! Loop index variable
!   
!  Calculations
!
   expon_approx: DO idx = 1, pc_group
!
      x_tilde(idx) = exp(-1.0_dp*lambda_l(idx)*delta_t)
      y_tilde(idx) = lambda_l(idx)*exp(-1.0_dp*lambda_l(idx)*delta_t)
!
      z_tilde(idx) = lambda_l(idx) * exp(-1.0_dp*(omega_p + lambda_l(idx)) * delta_t)* &
                   ((exp((omega_p+lambda_l(idx))*delta_t) - 1.0_dp) / (omega_p + lambda_l(idx)) - 1.0_dp)
!
   END DO expon_approx
END SUBROUTINE pc_expon_approx
!=======================================================================================================================
!
!  ************
!  * get_hats *
!  ************
!
!  Purpose: Subroutine to calculate x_hat, y_hat, and z_hat for n(i+1) equations
!           in the analytical precursor group integration.
!=======================================================================================================================
SUBROUTINE get_hats(&
   gen_time,        & ! IN (s)
   lambda_l,        & ! IN (1/s)
   c_bos,           & ! IN (1/cm^3)
   beta_l,          & ! IN (-)
   x_tilde,         & ! IN (-)
   y_tilde,         & ! IN (-)
   z_tilde,         & ! IN (-)
   x_hat,           & ! OUT (-)
   y_hat,           & ! OUT (-)
   z_hat,           & ! OUT (-)
   c_hat)             ! OUT (-)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
! Arguments
!
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Precursor decay constants  (1/s)
   REAL(dp), INTENT(IN) :: gen_time    ! Mean neutron generation time 
   REAL(dp), INTENT(IN) :: c_bos(:)    ! Beginning of step precursors' concentration at t(i) (1/cm^3)
   REAL(dp), INTENT(IN) :: beta_l(:)   ! Group-wise core delayed neutron fractions  (no unit)       
   REAL(dp), INTENT(IN) :: x_tilde(:)  ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(IN) :: y_tilde(:)  ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(IN) :: z_tilde(:)  ! Factor in EOS delayed n. pre. density equation	  
!
   REAL(dp), INTENT(OUT) :: x_hat ! Factor in EOS neutron density equation
   REAL(dp), INTENT(OUT) :: y_hat ! Factor in EOS neutron density equation
   REAL(dp), INTENT(OUT) :: z_hat ! Factor in EOS neutron density equation
   REAL(dp), INTENT(OUT) :: c_hat ! Factor in EOS neutron density equation
!
!  Locals
!
   INTEGER :: idx ! Loop index variable
!
!  Calculations
!
!  Set these to 0. In case of any previous calls in the iteration, wrong sums may be obtained otherwise.
   z_hat = 0.0_dp
   y_hat = 0.0_dp
   x_hat = 0.0_dp
   c_hat = 0.0_dp
!   
!  Calculate z_hat
!  ---------------
   zhat: DO idx = 1, pc_group
      z_hat = z_hat + z_tilde(idx) * beta_l(idx) / gen_time
   ENDDO zhat
!
!  Calculate y_hat
!  ---------------
   yhat: DO idx = 1, pc_group
      y_hat = y_hat + y_tilde(idx) * beta_l(idx) / gen_time
   ENDDO yhat
!
!  Calculate x_hat
!  ---------------
   xhat: DO idx = 1, pc_group
      x_hat = x_hat + lambda_l(idx) * x_tilde(idx) * c_bos(idx)
   ENDDO xhat
!  
!   Calculate c_hat
!   ---------------
   chat: DO idx = 1, pc_group
      c_hat = c_hat + lambda_l(idx) * c_bos(idx)
   ENDDO chat
!
END SUBROUTINE get_hats
!=======================================================================================================================
!
!  ****************
!  * update_n_eos *
!  ****************
!
!  Purpose: Subroutine to calculate n_eos.
!
!=======================================================================================================================
SUBROUTINE update_n_eos(&
   gen_time,            & ! IN (s)
   lambda_l,            & ! IN (1/s)
   delta_t,             & ! IN (s)
   rho_eos,             & ! IN (-)
   rho_bos,             & ! IN (-)
   beta_l,              & ! IN (-)
   beta,                & ! IN (-)
   x_tilde,             & ! IN (-)
   y_tilde,             & ! IN (-)
   z_tilde,             & ! IN (-)
   omega_n,             & ! IN (-)
   theta,               & ! IN (-)
   c_bos,               & ! IN (-)
   n_bos,               & ! IN (1/cm^3)
   n_eos)                 ! OUT (1/cm^3)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: gen_time    ! Mean neutron generation time (s)
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Precursor decay constants  (1/s)
   REAL(dp), INTENT(IN) :: delta_t     ! Time step size (s)
   REAL(dp), INTENT(IN) :: rho_eos     ! Reactivity at the end of step  
   REAL(dp), INTENT(IN) :: rho_bos     ! Reactivity at the beginning of step
   REAL(dp), INTENT(IN) :: beta        ! Delayed neutron fraction
   REAL(dp), INTENT(IN) :: beta_l(:)   ! Group-wise delayed neutron fractions        
   REAL(dp), INTENT(IN) :: x_tilde(:)  ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(IN) :: y_tilde(:)  ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(IN) :: z_tilde(:)  ! Factor in EOS delayed n. pre. density equation
   REAL(dp), INTENT(IN) :: n_bos       ! Beginning of step neutron density at t(i) (1/cm^3)
   REAL(dp), INTENT(IN) :: omega_n     ! Neutron density transformation frequency (no unit)
   REAL(dp), INTENT(IN) :: theta       ! Needed for theta method
   REAL(dp), INTENT(IN) :: c_bos(:)    ! Beginning of cycle precursor equations
!
   REAL(dp), INTENT(OUT) :: n_eos   ! End of step core average neutron density at t(i) (1/cm^3)
!
!  Locals
!
   REAL(dp) :: A        ! Factor in EOS neutron density equation
   REAL(dp) :: B        ! Factor in EOS neutron density equation
   REAL(dp) :: Q
   REAL(dp) :: freq_exp ! Factor a1 and a2 equation
   REAL(dp) :: x_hat    ! Factor in EOS neutron density equation
   REAL(dp) :: y_hat    ! Factor in EOS neutron density equation
   REAL(dp) :: z_hat    ! Factor in EOS neutron density equation
   REAL(dp) :: c_hat    ! Factor in EOS neutron density equation
   REAL(dp) :: test 
!      
!  Calculations
! 
   CALL get_hats(&
      gen_time = gen_time, & ! IN (1/s)
      lambda_l = lambda_l, & ! IN (s)
      c_bos    = c_bos   , & ! IN (1/cm^3)
      beta_l   = beta_l  , & ! IN (-)
      x_tilde  = x_tilde , & ! IN (-)
      y_tilde  = y_tilde , & ! IN (-)
      z_tilde  = z_tilde , & ! IN (-)
      x_hat    = x_hat   , & ! OUT (-)
      y_hat    = y_hat   , & ! OUT (-)
      z_hat    = z_hat   , & ! OUT (-)
      c_hat    = c_hat)      ! OUT (-)
!
!  Calculate freq_exp, a2 and a1
!  -----------------------------
!   For testing the frequencies.
!   c_hat = 0.0_dp
!   z_hat = 0.0_dp
!   x_hat = 0.0_dp
!   y_hat = 0.0_dp
   freq_exp = exp(-1.0_dp * omega_n * delta_t)
   test = (rho_bos - beta) / gen_time - omega_n
   B = freq_exp *(1.0_dp - delta_t * theta * ((rho_eos - beta) / gen_time + z_hat - omega_n ))
   A = 1.0_dp+delta_t * (1.0_dp - theta) * ((rho_bos - beta) / gen_time - omega_n) + theta * freq_exp * delta_t * y_hat
   Q =  delta_t * ( (1.0_dp -theta) * c_hat + freq_exp * theta * x_hat )
!
!  Calculate n(t_i+1)
!  ------------------
   n_eos = n_bos * A / B + delta_t * (1.0_dp -theta) * c_hat / B + freq_exp * delta_t * theta * x_hat / B
   109 FORMAT (ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6, ES14.6)
   WRITE(*,109) n_bos, n_eos, omega_n, freq_exp, test, A, B, Q, A/B, Q/B, z_hat
!   108 FORMAT (ES14.6)
!   WRITE(*,108) y_hat
!  
END SUBROUTINE update_n_eos
!=======================================================================================================================
!
!  ****************
!  * update_c_eos *
!  ****************
!
!  Purpose: Subroutine to calculate c_eos.
!
!=======================================================================================================================
SUBROUTINE update_c_eos(&
   gen_time,            & ! IN (s)
   lambda_l,            & ! IN (1/s)
   beta_l,              & ! IN (-)
   x_tilde,             & ! IN (-)
   y_tilde,             & ! IN (-)
   z_tilde,             & ! IN (-)
   n_bos,               & ! IN (1/cm^3)
   n_eos,               & ! IN (1/cm^3)
   c_bos,               & ! IN (1/cm^3)
   c_eos)                 ! OUT (1/cm^3)
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
!
!  Arguments
!
   REAL(dp), INTENT(IN) :: gen_time    ! Mean neutron generation time          (s)
   REAL(dp), INTENT(IN) :: lambda_l(:) ! Table to store precursor decay constants (1/s)
   REAL(dp), INTENT(IN) :: beta_l(:)   ! Delayed neutron fractions       
   REAL(dp), INTENT(IN) :: x_tilde(:)  ! EOS delayed n. pre. density equation factor
   REAL(dp), INTENT(IN) :: y_tilde(:)  ! EOS delayed n. pre. density equation factor
   REAL(dp), INTENT(IN) :: z_tilde(:)  ! EOS delayed n. pre. density equation factor 
   REAL(dp), INTENT(IN) :: n_bos       ! Beginning of step core average neutron density at t(i) (1/cm^3)
   REAL(dp), INTENT(IN) :: n_eos       ! End of step core average neutron density
   REAL(dp), INTENT(IN) :: c_bos(:)    ! Beginning of step core average precursor concentrations
!  
   REAL(dp), INTENT(OUT) :: c_eos(:)   ! End of step core average precursor concentration
!
!  Locals
!
   INTEGER  :: idx              ! Loop index variable
   REAL(dp) :: Gl_bar(pc_group) ! Factor in the c_eos equation
!      
! Calculations
!
!  For each precursor group calculate Gl and c_eos
!  -----------------------------------------------
   c_eos_update: DO idx = 1, pc_group
      Gl_bar(idx) = beta_l(idx) / gen_time
      c_eos(idx) = x_tilde(idx) * c_bos(idx) + y_tilde(idx) / lambda_l(idx) * Gl_bar(idx) * n_bos &
                   + z_tilde(idx) / lambda_l(idx) * Gl_bar(idx) * n_eos
   ENDDO c_eos_update
!  
END SUBROUTINE update_c_eos
!=======================================================================================================================
!
!  *******************
!  * pk_methods_init *
!  *******************
!
!  Purpose: Subroutine to initilize the pk_methods module.
!
!=======================================================================================================================
SUBROUTINE pk_methods_init()
!=======================================================================================================================
! Record of revisions:
!        Date       Programmer               Description of changes
!        -----------------------------------------------------------
!        24/05/19   K. Luszczek              Original code
!
!=======================================================================================================================
!  Check if only one pass through the main loop is needed.
!  Used if either neutron density frequency transformation is to be used 
!  or if the exponential fission source approximation is to be used in the precursor equation
!  If yes, then set max_n_eos_iter to 1.
!  ------------------------------------------------------------------------------------------
   omega = 0.0_dp ! Sets initial value of omega to zero.
!
   IF (n_solution == N_NO_FREQUENCY_TRANS) THEN
      IF (pre_approx /= C_EXPONENT_APPROX) THEN
         CALL set_data("max_n_eos_iter", 1)
      ENDIF
   ENDIF
!
!  If convergence on omega is requested. Omega HAS to be updated between iterations.
!
   IF (converg_opt == CONV_W) THEN
      CALL set_data("omega_option",1)
   ENDIF
!
END SUBROUTINE pk_methods_init
!=======================================================================================================================
END MODULE pk_methods
!=======================================================================================================================