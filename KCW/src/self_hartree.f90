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
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE buffers,              ONLY : get_buffer
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : nspin_mag
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
  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr, nrho)
  COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
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
   !! sh = sh + 0.5D0 * sum (CONJG(rhog (:,1)) * vh_rhog(:) )*weight(iq)*omega
    zpom = (0.0_dp, 0.0_dp)
    !$acc parallel loop reduction(+: zpom) present(rhog, vh_rhog)
    DO ii = 1, ngms
       zpom = zpom + CONJG(rhog (ii,1)) * vh_rhog(ii) 
    END DO   
    sh = sh + 0.5D0 * zpom*weight(iq)*omega
    !
    ! 
  ENDDO ! qpoints
  !$acc exit data delete(rhor, rhog , delta_vg, vh_rhog, delta_vg_ )
  !
  DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
  !
  CALL mp_sum (sh, intra_bgrp_comm)
 !
END SUBROUTINE self_hartree

