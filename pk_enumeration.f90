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
   INTEGER, PARAMETER, PUBLIC :: N_FULLY_IMPLICIT  = 0 ! Use fully implicit backward Euler for 
!                                                        the neutron density equation
   INTEGER, PARAMETER, PUBLIC :: N_FREQUENCY_TRANS = 1 ! Use frequency transformation and consecutive theta integration
!                                                        for the neutron density equation
!
!  Enumeration for plotting
!
   INTEGER, PARAMETER, PUBLIC :: PLOT_RHO = 2 ! Corresponds to the column number in output_data.dat.
   INTEGER, PARAMETER, PUBLIC :: PLOT_N   = 3 ! Corresponds to the column number in output_data.dat.
   INTEGER, PARAMETER, PUBLIC :: PLOT_C   = 4 ! Corresponds to the column number in output_data.dat.
!
!=======================================================================================================================
END MODULE pk_enumeration
!=======================================================================================================================
