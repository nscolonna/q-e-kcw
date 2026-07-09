!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
SUBROUTINE self_hartree (iwann, sh)
  !----------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : num_wann, nqstot, iurho_wann, nrho
  USE fft_base,             ONLY : dffts
  USE cell_base,            ONLY : omega
  USE gvecs,                ONLY : ngms
  USE gvect,                ONLY : gstart
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE buffers,              ONLY : get_buffer
  USE noncollin_module,     ONLY : nspin_mag
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! The scalar contribution to the hamiltonian 
  INTEGER, INTENT(IN) :: iwann
  !
  ! Couters for the q point, wannier index. record length for the wannier density
  INTEGER :: iq, lrrho, ii
  !
  ! The periodic part of the wannier orbital density
  !COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr, nrho)
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:,:), rhor(:,:)
  !COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
  COMPLEX(DP), ALLOCATABLE :: delta_vr(:,:), delta_vr_(:,:)
  !
  ! The self Hartree
  COMPLEX(DP) :: sh, zpom
  !
  ! Auxiliary variables 
  COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
  !
  ! The weight of each q point
  REAL(DP) :: weight(nqstot)
  !
  ALLOCATE ( rhog (ngms,nrho) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
  ALLOCATE ( rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr, nrho) )
  ALLOCATE ( delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag) )
  !
  !$acc enter data create(rhor, rhog, vh_rhog, delta_vr, delta_vr_, delta_vg, delta_vg_)
  DO iq = 1, nqstot
    !
    lrrho=num_wann*dffts%nnr*nrho
    CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    !! Retrive the rho_wann_q(r) from buffer in REAL space
    !
    weight(iq) = 1.D0/nqstot ! No SYMM 
    !
   !$acc kernels present(rhog, delta_vg, delta_vg_, vh_rhog, rhor)
    rhog(:,:)       = (0.0_dp,0.0_dp)
    delta_vg(:,:)   = (0.0_dp,0.0_dp)
    vh_rhog(:)      = (0.0_dp,0.0_dp)
    rhor(:,:)       = (0.0_dp,0.0_dp)
    !$acc end kernels
    !
    rhor(:,:) = rhowann(:,iwann,:) 
    !$acc update device(rhor)
    !! The periodic part of the orbital desity in real space
    !
    CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
    !! The periodic part of the perturbation DeltaV_q(G)
    ! 
    IF (gamma_only) THEN 
      sh = sh + DBLE ( sum (CONJG(rhog (:,1)) * vh_rhog(:) ) ) * weight(iq)*omega
      IF (gstart == 2) sh = sh - 0.5D0 * DBLE( CONJG(rhog (1,1)) * vh_rhog(1) ) *weight(iq)*omega
    ELSE
      sh = sh + 0.5D0 * sum (CONJG(rhog (:,1)) * vh_rhog(:) )*weight(iq)*omega
    ENDIF
    !
    ! 
  ENDDO ! qpoints
  !$acc exit data delete(rhor, rhog , delta_vg, vh_rhog, delta_vg_ )
  !$acc exit data delete(delta_vr, delta_vr_)
  !
  DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
  DEALLOCATE ( rhowann, rhor )
  DEALLOCATE ( delta_vr, delta_vr_ )
  !
  CALL mp_sum (sh, intra_bgrp_comm)
 !
END SUBROUTINE self_hartree

