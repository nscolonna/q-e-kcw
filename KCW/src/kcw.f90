!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
PROGRAM kcw
  !-----------------------------------------------------------------
  !
  !!  This is the main driver of the kcw.x code. It computes the 
  !!  screaning parameters alpha for KOOPMANS. It reads the ouput of 
  !!  a PWscf calculation and of wann2kcw and computes
  !!  the orbital dependent screening coefficients as described in 
  !!  N. Colonna et al. JCTC 14, 2549 (2018) 
  !!  https://pubs.acs.org/doi/10.1021/acs.jctc.7b01116
  !!  The KC Hamiltonian is also computed, possibly iterpolated
  !!  and diagonalized.
  !!
  !!  Code written by Nicola Colonna (EPFL April 2019) 
  !
  USE mp_global,         ONLY : mp_startup,mp_global_end 
  USE environment,       ONLY : environment_start, environment_end
  USE control_kcw,       ONLY : calculation
  USE mp_global,         ONLY : mp_startup
  USE check_stop,        ONLY : check_stop_init
  USE coulomb,           ONLY : setup_coulomb
  USE io_global,         ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code='KCW'
  !
  CALL header ()
  !
  ! 1) Initialize MPI, clocks, print initial messages
  CALL mp_startup ( )
  !
  CALL environment_start ( code )
  !
  !
  CALL check_stop_init ()
  !
  ! 2) Read the input file and the PW outputs
  CALL kcw_readin( ) 
  !
  IF (calculation == 'cc') call setup_coulomb()
  IF (calculation == 'wann2kcw') CALL wann2kcw ( )
  IF (calculation == 'screen')   CALL kcw_screen ( )
  IF (calculation == 'ham' )     CALL kcw_ham ()
  !
  CALL print_clock_kcw ( )
  !
  ! 5) Clean and Close 
  CALL mp_global_end()
  CALL environment_end( code )
  !
END PROGRAM kcw


SUBROUTINE header 
  !
  USE io_global,         ONLY : stdout
  IMPLICIT NONE

  WRITE( stdout, '(/5x,"=--------------------------------------------------------------------------------=")')
  WRITE(stdout,*) "                     :::    :::           ::::::::         :::       ::: "
  WRITE(stdout,*) "                    :+:   :+:           :+:    :+:        :+:       :+:  "
  WRITE(stdout,*) "                   +:+  +:+            +:+               +:+       +:+   "
  WRITE(stdout,*) "                  +#++:++             +#+               +#+  +:+  +#+    "
  WRITE(stdout,*) "                 +#+  +#+            +#+               +#+ +#+#+ +#+     "
  WRITE(stdout,*) "                #+#   #+#           #+#    #+#         #+#+# #+#+#       " 
  WRITE(stdout,*) "               ###    ###           ########           ###   ###         "
  WRITE( stdout, '(/5x,"  Koopmans functional implementation based on DFPT; please cite this program as")')
  WRITE( stdout, '(/5x,"     N.Colonna, R. De Gannaro, E. Linscott, and N. Marzari, arXiv:2202.08155   ")')
  WRITE( stdout, '( 5x,"=--------------------------------------------------------------------------------=")')
  !
END SUBROUTINE header
