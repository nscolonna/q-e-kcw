!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE ks_hamiltonian (evc, ik, h_dim, eigvl_out)
  !---------------------------------------------------------------------
  !
  !! This routine compute and diagonalize the KS Hamiltonian
  !! Non-collinear case is NOT implemented!
  !! OBSOLETE?
  !
  !! ik is the LOCAL (pool) k-point index. When check_ks is on, the caller
  !! passes eigvl_out to retrieve the "WANN" (KS-in-Wannier-gauge) eigenvalues
  !! computed here instead of having this routine print them directly: with
  !! pools active this is called once per LOCAL k-point, so a WRITE here would
  !! only ever reach the log for the k-points owned by ionode's own pool. The
  !! caller is expected to gather eigvl_out (and the corresponding et(:,ik))
  !! across pools and print the full, ordered table once - see kcw_setup_ham.f90
  !! and rotate_ks.f90.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : npwx, npw, et
  USE uspp,                 ONLY : nkb
  USE becmod,               ONLY : becp, allocate_bec_type_acc, deallocate_bec_type_acc
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
  USE klist,                ONLY : init_igk, xk, nkstot
  USE mp,                   ONLY : mp_sum
  USE control_kcw,          ONLY : Hamlt, calculation, spin_component, check_ks
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_lsda, nspin_gga, nspin_mag

  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)    :: ik, h_dim
  !
  COMPLEX(DP), INTENT(IN) :: evc(npwx*npol,h_dim)
  !
  REAL(DP), INTENT(OUT), OPTIONAL :: eigvl_out(h_dim)
  ! The "WANN" (KS-in-Wannier-gauge) eigenvalues, returned to the caller when check_ks is on
  !
  !! COMPLEX(DP) :: hpsi(npwx*npol,h_dim), ham(h_dim,h_dim), hij, eigvc(npwx*npol,h_dim)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), ham(:,:), eigvc(:,:)
  COMPLEX(DP) :: hij
  !
  !! REAL(DP) :: eigvl(h_dim), check
  REAL(DP) :: check
  REAL(DP), ALLOCATABLE :: eigvl(:)
  !
  INTEGER :: iband, jband, ig, ik_eff
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  !
  CALL allocate_bec_type_acc ( nkb, h_dim, becp, intra_bgrp_comm )
  !
  CALL init_igk ( npwx, ngm, g, gcutw )
  !
  call g2_kin( ik )
  ALLOCATE(hpsi(npwx*npol,h_dim), ham(h_dim,h_dim),eigvc(npwx*npol,h_dim), eigvl(h_dim))
  !$acc enter data create(hpsi) copyin(evc)
  !$acc kernels present(hpsi)
  hpsi(:,:) = (0.D0, 0.D0)
  !$acc end kernels
  !
  CALL h_psi( npwx, npw, h_dim, evc, hpsi )
  !
  ! ##### Build up the KI Hamiltonian 
  !
  !$acc data copyout(ham) present(hpsi, evc)
  !$acc kernels
  ham(:,:)= (0.D0,0D0)
  !$acc end kernels
  !
  
  DO iband = 1, h_dim
     ! 
     !$acc parallel loop private(hij)
     DO jband = iband, h_dim
        !
        hij = (0.0_dp, 0.0_dp)
        !$acc loop reduction(+:hij)
        DO ig = 1, npw*npol
           hij = hij + CONJG(evc(ig,iband)) * hpsi(ig,jband)
        ENDDO
        !! CALL mp_sum (hij, intra_bgrp_comm)
        !
        ham(iband,jband) = hij
        ham(jband,iband) = CONJG(ham(iband,jband))
        !
        !IF (iband==jband) WRITE(*,*) iband,jband, REAL(hij)*rytoev
     ENDDO
     !
  ENDDO
  !$acc end data
  !$acc exit data  delete(evc, hpsi) 
!! sum ouside the loops
   CALL mp_sum (ham, intra_bgrp_comm)
   
  !
  ! Store the hamiltonian in the Wannier Gauge
  !
  IF (calculation == 'ham') then
    ik_eff = global_kpoint_index (nkstot, ik) - (spin_component -1)*nkstot/nspin
    !WRITE(*,*) ik, ik_eff
    Hamlt(ik_eff,1:h_dim,1:h_dim) = ham(1:h_dim,1:h_dim)
  ENDIF
  !
  ! Check the eigenvalue are consistent with the PWSCF calculation
  IF (check_ks) THEN
    CALL cdiagh( h_dim, ham, h_dim, eigvl, eigvc )
  ENDIF
  !
  check = 0.D0
  DO iband = 1, h_dim
    check = check + (eigvl(iband)-et(iband,ik))/h_dim
  ENDDO 
  !
  IF ( check_ks .AND. PRESENT(eigvl_out) ) eigvl_out(1:h_dim) = eigvl(1:h_dim)
  ! The caller prints the WANN/PWSCF report (see note above)
  !
  CALL deallocate_bec_type_acc (becp)

  DEALLOCATE(hpsi, ham, eigvc, eigvl)
  !
END SUBROUTINE ks_hamiltonian
