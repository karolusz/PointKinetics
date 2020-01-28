!=======================================================================================================================
!
!  ******************
!  * pk_enumeration *
!  ******************
!
!  Purpose:
!  --------
!
!  This module contains enumerations used in the point kinetics module
!
!=======================================================================================================================
MODULE pk_enumeration
!=======================================================================================================================
! Record of revisions:
!       Date       Programmer               Description of changes
!       -----------------------------------------------------------
!       27/09/18   K. Luszczek              Original code
!
!=======================================================================================================================
   USE pk_kinds     , ONLY: dp
!
   IMPLICIT NONE
!
!  Enumeration for a method of precursor equation solution
!
   INTEGER, PARAMETER, PUBLIC :: C_FULLY_IMPLICIT  = 0 ! Use fully implicit backward Euler in precursor eq.
   INTEGER, PARAMETER, PUBLIC :: C_CONSTANT_APPROX = 1 ! Use constant fission source approximation in precursor eq.
   INTEGER, PARAMETER, PUBLIC :: C_LINEAR_APPROX   = 2 ! Use linear fission source approximation in precursor eq.
   INTEGER, PARAMETER, PUBLIC :: C_EXPONENT_APPROX = 3 ! Use exponential fission source approximation in precursor eq.
! 
!  Enumeration for a method of neutron density equation solution
!
   INTEGER, PARAMETER, PUBLIC :: N_NO_FREQUENCY_TRANS = 0 ! Do not use the frequency transfomration for the n. dens. eq.
   INTEGER, PARAMETER, PUBLIC :: N_FREQUENCY_TRANS    = 1 ! Use frequency transformation for the neutron density equation
!
!  Enumeration for plotting
!
   INTEGER, PARAMETER, PUBLIC :: PLOT_RHO = 2 ! Corresponds to the column number in output_data.dat.
   INTEGER, PARAMETER, PUBLIC :: PLOT_N   = 3 ! Corresponds to the column number in output_data.dat.
   INTEGER, PARAMETER, PUBLIC :: PLOT_C   = 4 ! Corresponds to the column number in output_data.dat.
!
!  Enumeration for updating the omega during iterations
!
   INTEGER, PARAMETER, PUBLIC :: OMEGA_NO  = 0 ! Omega is not updated between iterations. Only between steps.
   INTEGER, PARAMETER, PUBLIC :: OMEGA_YES = 1 ! Omega is updated between iterations and between steps.
!
!  Enumeration for convergence option (neutron denisty/ omega frequency)
!
   INTEGER, PARAMETER, PUBLIC :: CONV_N = 0 ! The inner loop will converge on n.
   INTEGER, PARAMETER, PUBLIC :: CONV_W = 1 ! The inner lopp will converge on omega.
!
!  Enumeration for the time step control option
!
   INTEGER, PARAMETER, PUBLIC :: DT_ADAPT = 0 ! Calculate adaptive time steps.
   INTEGER, PARAMETER, PUBLIC :: DT_CONST = 1 ! Use constant time step from user input.
!
!=======================================================================================================================
END MODULE pk_enumeration
!=======================================================================================================================
