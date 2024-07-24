!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!adapted from PH code
!-----------------------------------------------------------------------
SUBROUTINE kcw_psymdvscf (dvtosym)!(nper, irr, dvtosym)
  !-----------------------------------------------------------------------
  !! p-symmetrize the potential.
  !
  USE kinds,     ONLY : DP
  USE noncollin_module, ONLY : nspin_mag
  USE mp_bands,  ONLY : me_bgrp
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym

  USE lr_symm_base, ONLY : nsymq, minus_q
  USE control_kcw,  ONLY : irr_bz
  !
  IMPLICIT NONE
  !
  !INTEGER :: nper
  !! the number of perturbations
  !INTEGER :: irr
  !! the representation under consideration
  COMPLEX(DP) :: dvtosym(dfftp%nnr,nspin_mag)!,nper)
  !! the potential to symmetrize
  !
  ! ... local variable
  !
#if defined (__MPI)
  !
  INTEGER :: i, is, iper, ir3, ioff, ioff_tg, nxyp

  COMPLEX(DP), ALLOCATABLE :: ddvtosym (:,:)!,:)
  ! the potential to symm

  if ( .not. irr_bz ) RETURN
  IF (nsymq.EQ.1.AND. (.NOT.minus_q) ) RETURN
  CALL start_clock ('psymdvscf')

  ALLOCATE (ddvtosym ( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin_mag) )!, nper))

  !DO iper = 1, nper
     DO is = 1, nspin_mag
        CALL cgather_sym (dfftp, dvtosym (:, is), ddvtosym(:, is) ) !, iper), ddvtosym (:, is) )!, iper) )
     ENDDO

  !ENDDO

  CALL kcw_symdvscf (ddvtosym) !(nper, irr, ddvtosym)

  nxyp = dfftp%nr1x * dfftp%my_nr2p
  !DO iper = 1, nper
     DO is = 1, nspin_mag
        DO ir3 = 1, dfftp%my_nr3p
           ioff    = dfftp%nr1x * dfftp%my_nr2p * (ir3-1)
           ioff_tg = dfftp%nr1x * dfftp%nr2x    * (dfftp%my_i0r3p+ir3-1) + dfftp%nr1x * dfftp%my_i0r2p 
           CALL zcopy (nxyp, ddvtosym (ioff_tg+1, is), 1, dvtosym (ioff+1, is), 1)
        END DO
     ENDDO
  !ENDDO
  DEALLOCATE (ddvtosym)

  CALL stop_clock ('psymdvscf')
#else
   if ( .not. irr_bz ) RETURN
  CALL kcw_symdvscf (dvtosym)!(nper, irr, dvtosym)
#endif

  RETURN

END SUBROUTINE kcw_psymdvscf
