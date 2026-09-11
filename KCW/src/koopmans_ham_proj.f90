! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!-----------------------------------------------------------------------
SUBROUTINE koopmans_ham_proj (delta)
  !---------------------------------------------------------------------
  ! Here the KI hamiltonian is written in terms of projectors on Wannnier 
  ! functions:
  ! \Delta_H_KI = \sum_n (1/2-P_n) * \Delta_n * |w_kn><w_kn| with 
  !   - P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> the occupation of the Wannier function n, 
  !   - \Delta_n=\alpha_n <W_n^2 | f_Hxc | W_n^2> the onsite KI correction 
  ! This correction can be applied perturbatively or not (default) on all the KS states
  ! available from the preceeding nscf calculation
  !
  USE io_global,             ONLY : stdout, ionode
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : xk, ngk, igk_k, nkstot, nks
  USE control_kcw,           ONLY : num_wann, evc0, Hamlt, kcw_iverbosity, &
                                    num_wann_occ, iuwfc_wann_allk, spin_component, nkstot_eff
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd, current_k
  USE units_lr,              ONLY : iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum, mp_max, mp_min
  USE mp_pools,              ONLY : inter_pool_comm
  USE lsda_mod,              ONLY : lsda, isk, current_spin
  USE uspp,                  ONLY : nkb, vkb
  USE uspp_init,             ONLY : init_us_2
  USE noncollin_module,      ONLY : npol
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !! The global index of a local (pool) k-point
  !
  ! ik_loc is the LOCAL (pool) k-point index: the KS orbitals (iuwfc), the number of
  ! PWs (ngk), the eigenvalues (et), igk_k and xk are all pool-local arrays, so they
  ! must be addressed with it. ik is the "effective" (1:nkstot_eff) index of the
  ! current spin channel, used for the pool-replicated quantities: Hamlt and the
  ! ALL-k Wannier-gauge buffer iuwfc_wann_allk.
  INTEGER :: ik, ibnd, ik_loc, k, eig_start, eig_win, eig_top, nvals
  !
  ! Per-k results, stashed at the effective k index and gathered across pools after
  ! the loop so that the full, correctly ordered table is printed once from ionode.
  REAL(DP), ALLOCATABLE :: eigvl_ks_all(:,:), eigvl_ki_all(:,:), eigvl_pert_all(:,:)
  REAL(DP), ALLOCATABLE :: xk_all(:,:)
  REAL(DP), ALLOCATABLE :: proj_spec(:,:,:)
  ! kcw_iverbosity>1 only: KI eigenvalues around the Fermi level as a function of the
  ! size of the diagonalized Hilbert space
  !
  ! The on-site KI correction \alpha_n*<W_0n^2|f_Hxc|W_0n^2> (computed in dH_ki_wann.f90)
  COMPLEX(DP), INTENT(IN) :: delta (num_wann)
  ! and occupations numbers of WFs: P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> 
  REAL(DP) :: occ_mat(num_wann)
  !
  COMPLEX(DP) :: overlap
  COMPLEX(DP) :: ham(nbnd,nbnd), eigvc(nbnd,nbnd), deltah(nbnd,nbnd)
  !
  REAL(DP) :: eigvl_wann_chk(nbnd)
  ! Throwaway receiver for the (mandatory) eigvl_out argument of ks_hamiltonian
  !
  ! The new eigenalues 
  REAL(DP), ALLOCATABLE :: eigvl_aux(:)
  COMPLEX(DP) , ALLOCATABLE :: eigvc_aux(:,:), ham_aux(:,:), evc_aux(:,:)
  !
  ! The new eigenalues
  REAL(DP) :: eigvl(nbnd)
  REAL(DP) :: eigvl_pert(nbnd)
  REAL(DP) :: eigvl_ks(nbnd)
  !
  INTEGER :: i
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  REAL(DP) :: ehomo_pert, elumo_pert
  INTEGER  :: lrwannfc
  REAL(DP), EXTERNAL :: get_clock
  !
  ALLOCATE (evc_aux(npwx*npol,nbnd))
  !
  WRITE( stdout, '(/,5X, "INFO: BUILD and DIAGONALIZE the KI HAMILTONIAN")')
  WRITE( stdout, '(  5X, "INFO: Projectors scheme")')
  !
  ! The occupation matrix
  ! P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> 
  CALL occupations(occ_mat)
  WRITE(stdout, 900) get_clock('KCW')
  ! 
#ifdef DEBUG
  WRITE(stdout,'(/,5X,"Screened on-site correction:")')
  WRITE(stdout,'(5X,2(F10.6, 2x), F10.6)') (delta(iwann), occ_mat(iwann), iwann=1, num_wann)
#endif
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  ehomo_pert=-1D+6
  elumo_pert=+1D+6
  !
  ALLOCATE ( eigvl_ks_all(nbnd, nkstot_eff), eigvl_ki_all(nbnd, nkstot_eff), &
             eigvl_pert_all(nbnd, nkstot_eff), xk_all(3, nkstot_eff) )
  eigvl_ks_all = 0.D0; eigvl_ki_all = 0.D0; eigvl_pert_all = 0.D0; xk_all = 0.D0
  !
  ! ... Size of the kcw_iverbosity>1 report (see the print block after the loop)
  eig_win   = 5
  eig_start = MAX(num_wann_occ-eig_win+1, 1)
  eig_top   = num_wann_occ + eig_win
  nvals     = MIN(nbnd, eig_top) - eig_start + 1
  IF (kcw_iverbosity .gt. 1 .AND. nvals .gt. 0) THEN
    ALLOCATE ( proj_spec(nvals, nbnd, nkstot_eff) )
    proj_spec = 0.D0
  ENDIF
  !
  ! ... Loop over the LOCAL (this pool's) k-points only: each pool can only read the
  ! KS orbitals of the k-points it owns from its own iuwfc buffer.
  !
  DO ik_loc = 1, nks
    !
    IF ( lsda .AND. isk(ik_loc) /= spin_component ) CYCLE
    !
    ik = global_kpoint_index (nkstot, ik_loc) - (spin_component-1)*nkstot_eff
    !
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_loc )
    npw = ngk(ik_loc)
    xk_all(:,ik) = xk(:,ik_loc)
    !
    ehomo_ks = MAX ( ehomo_ks, et(num_wann_occ  , ik_loc) )
    IF (nbnd > num_wann_occ) elumo_ks = MIN ( elumo_ks, et(num_wann_occ+1, ik_loc) )
    !
    !
    ! Build and diagonalize the projector-based KI Hamiltonian on the Hilbert space spanned
    ! by the KS states available from the preceeding nscf calculation
    ! KI contribution at k: deltah_ij = \sum_n [ (1/2-P_n)D_n <u_ki | w_kn><w_kn | u_kj>]
    !
    IF (.FALSE.) THEN
      ! In the canonicla KS basis the KS hamiltonian is already diagonal
      ! maening there is no need to re-build it. I keep this as is in case
      ! we want to use a diffferent basis
      current_k = ik_loc
      IF ( lsda ) current_spin = isk(ik_loc)
      IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik_loc), xk(1,ik_loc), vkb )
      ! FIXME: ned o modify ks_hamiltonian to pass the hammiltonian (and not use Hamlt of kcw_comm)
      ! NB: the 4th argument is the (mandatory) eigvl_out of ks_hamiltonian; it used to
      ! be called here with a LOGICAL .false. in that slot, which silently mismatched the
      ! dummy argument list (no explicit interface catches it). Unused here: check_ks
      ! is not active on this path, so ks_hamiltonian leaves it untouched.
      CALL ks_hamiltonian (evc, ik_loc, nbnd, eigvl_wann_chk)
      !
      ! The KS hamiltonian in the Wannier Gauge (just to check)
      ham(:,:) = Hamlt(ik,:,:) 
      CALL cdiagh( nbnd, ham, nbnd, eigvl, eigvc )
      !
      !WRITE( stdout, 9020 ) ( xk(i,ik_loc), i = 1, 3 )
      ! this shoud perfetcly matchs with the KS eigenvalues. If not, there is a problem
      WRITE( stdout, '(10x, "KS* ",8F11.4)' ) (eigvl(ibnd)*rytoev, ibnd=1,nbnd)
      !
    ELSE
      !
      ! The KS Hamiltonian in the KS basis
      ham(:,:)=CMPLX(0.D0, 0.D0, kind=DP)
      DO i = 1, nbnd 
        ham(i,i)    = et(i,ik_loc)
        eigvl_ks(i) = et(i,ik_loc)
      ENDDO
      !
    ENDIF
    !
    ! \Delta_H_KI_ij = <psi_i | \Delta_H_KI | psi_j> with
    ! \Delta_H_KI = \sum_n (1/2-P_n) * \Delta_n * |w_kn><w_kn|
    CALL dki_hamiltonian (evc, ik, nbnd, occ_mat, delta, deltah) 
    !
    !
    ! Because we have defined a uniq KI Hamiltonian, we can do a perturbative approach
    ! i.e. we keep only the diagonal part of the KI Hamiltoniana
    DO i = 1, nbnd
      eigvl_pert(i) = et(i,ik_loc) + DBLE(deltah(i,i))
    ENDDO
    ehomo_pert = MAX ( ehomo_pert, eigvl_pert(num_wann_occ ) )
    IF (nbnd > num_wann_occ) elumo_pert = MIN ( elumo_pert, eigvl_pert(num_wann_occ+1 ) )
    !
    ! Add the KI contribution to the KS Hamiltonian
    ham(:,:) = ham(:,:) + deltah(:,:) 
    ! And Diagonalize it 
    CALL cdiagh( nbnd, ham, nbnd, eigvl, eigvc )
    !
    IF ( ALLOCATED(proj_spec) ) THEN
      !
      ! Stash the spectrum; printed after the loop (see below)
      DO k = eig_start, nbnd
         !
         ALLOCATE (ham_aux(k,k), eigvl_aux(k), eigvc_aux(k,k))
         !ham_aux(1:k,1:k) = ham(1:k,1:k)
         ham_aux(1:k,1:k) = ham(1:k,1:k)
         !
         CALL cdiagh( k, ham_aux, k, eigvl_aux, eigvc_aux )
         !
         !neig_max = MIN (20, nbnd)
         proj_spec(1:MIN(k,eig_top)-eig_start+1, k, ik) = eigvl_aux(eig_start:MIN(k,eig_top))
         !
         DEALLOCATE (ham_aux)
         DEALLOCATE (eigvl_aux, eigvc_aux)
         !
      ENDDO
    ENDIF
    !
    !Overwrite et and evc
    et(1:nbnd, ik_loc) = eigvl(1:nbnd)
    ! MB
    ! This is different wrt koopmans_ham.f90:
    ! (1) the first dimension (row of A/evc) = npwx, not npw;
    ! (2) cannot use the same input matrix as output; need evc_aux
    CALL ZGEMM( 'N','N', npwx*npol, nbnd, nbnd, ONE, evc, npwx*npol, eigvc, nbnd, &
    ZERO, evc_aux, npwx*npol )
    evc(:,:) = evc_aux(:,:)
    CALL save_buffer ( evc, nwordwfc, iuwfc, ik_loc )
    !
    ehomo = MAX ( ehomo, eigvl(num_wann_occ ) )
    IF (nbnd > num_wann_occ) elumo = MIN ( elumo, eigvl(num_wann_occ+1 ) )
    !
    ! Stash the eigenvalues; the per-k report is printed after the loop (see below)
    eigvl_ks_all(:,ik)   = eigvl_ks(:)
    eigvl_ki_all(:,ik)   = eigvl(:)
    eigvl_pert_all(:,ik) = eigvl_pert(:)
    !
    WRITE(stdout, 901) get_clock('KCW')
    !
  ENDDO
  !
  ! ... Gather across pools: each pool has filled only the columns of the k-points it
  ! owns (the arrays were zeroed above and each effective k index is owned by exactly
  ! one pool), so a sum reconstructs the full table on every process.
  !
  CALL mp_sum ( eigvl_ks_all,   inter_pool_comm )
  CALL mp_sum ( eigvl_ki_all,   inter_pool_comm )
  CALL mp_sum ( eigvl_pert_all, inter_pool_comm )
  CALL mp_sum ( xk_all,         inter_pool_comm )
  IF ( ALLOCATED(proj_spec) ) CALL mp_sum ( proj_spec, inter_pool_comm )
  CALL mp_max ( ehomo_ks,   inter_pool_comm )
  CALL mp_max ( ehomo,      inter_pool_comm )
  CALL mp_max ( ehomo_pert, inter_pool_comm )
  CALL mp_min ( elumo_ks,   inter_pool_comm )
  CALL mp_min ( elumo,      inter_pool_comm )
  CALL mp_min ( elumo_pert, inter_pool_comm )
  !
  ! ... The per-k report, now that every process holds the full gathered table:
  ! print once, in k-point order, from ionode only.
  !
  IF ( ionode ) THEN
    DO ik = 1, nkstot_eff
      WRITE( stdout, 9020 ) ( xk_all(i,ik), i = 1, 3 )
      IF ( ALLOCATED(proj_spec) ) THEN
        WRITE(stdout,'( 12x, "KI Eigenvalues around Fermi as a function of the Hilbert space")')
        WRITE(stdout,'( 12x, "KI Eigenvaluse from", I5, " to", I5, " (i.e. +/-", I3, " around VBM),/")') &
                eig_start, num_wann_occ+eig_win, eig_win
        WRITE(stdout,'( 12x, " dim   eig1      eig2      ...")')
        DO k = eig_start, nbnd
          WRITE(stdout,'(12x, I3, 10F10.4)') k, &
                proj_spec(1:MIN(k,eig_top)-eig_start+1, k, ik)*rytoev
        ENDDO
        WRITE(stdout,*)
      ENDIF
      WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (eigvl_ks_all(ibnd,ik)*rytoev, ibnd=1,nbnd)
      WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (eigvl_ki_all(ibnd,ik)*rytoev, ibnd=1,nbnd)
      WRITE( stdout, '(10x, "pKI ",8F11.4)' ) (eigvl_pert_all(ibnd,ik)*rytoev, ibnd=1,nbnd)
    ENDDO
  ENDIF
  !
  IF ( elumo < 1d+6) THEN
    WRITE( stdout, 9042 ) ehomo_ks*rytoev, elumo_ks*rytoev
    WRITE( stdout, 9044 ) ehomo*rytoev, elumo*rytoev
    WRITE( stdout, 9046 ) ehomo_pert*rytoev, elumo_pert*rytoev
  ELSE
    WRITE( stdout, 9043 ) ehomo_ks*rytoev
    WRITE( stdout, 9045 ) ehomo*rytoev
    WRITE( stdout, 9047 ) ehomo_pert*rytoev
  END IF
  !
9043 FORMAT(/,8x, 'KS  highest occupied level (ev): ',F10.4 )
9042 FORMAT(/,8x, 'KS  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9045 FORMAT(  8x, 'KI  highest occupied level (ev): ',F10.4 )
9044 FORMAT(  8x, 'KI  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9047 FORMAT(  8x, 'pKI highest occupied level (ev): ',F10.4 )
9046 FORMAT(  8x, 'pKI highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
900 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
901 FORMAT('          total cpu time spent up to now is ',F10.1,' secs' )
  !
  DEALLOCATE (evc_aux)
  DEALLOCATE (eigvl_ks_all, eigvl_ki_all, eigvl_pert_all, xk_all)
  IF ( ALLOCATED(proj_spec) ) DEALLOCATE (proj_spec)
  !
  RETURN
  !
  CONTAINS
  !
  ! !----------------------------------------------------------------
  SUBROUTINE dki_hamiltonian (evc, ik, h_dim, occ_mat, delta, deltah)
    !----------------------------------------------------------------
    !
    USE buffers,               ONLY : get_buffer
    USE control_kcw,           ONLY : num_wann
    USE wvfct,                 ONLY : npwx
    USE control_flags,         ONLY : gamma_only
    USE gvect,                 ONLY : gstart
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: h_dim
    COMPLEX(DP) :: delta (num_wann)
    COMPLEX(DP), INTENT(IN) :: evc(npwx*npol,h_dim)
    INTEGER, INTENT(IN) :: ik
    REAL(DP), INTENT(IN) :: occ_mat(num_wann)
    COMPLEX(DP), INTENT(OUT) :: deltah(h_dim,h_dim)
    !
    INTEGER :: lrwannfc, ib, jb, iwann
    COMPLEX(DP) :: overlap_in, overlap_nj, overlap
    !
    lrwannfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann_allk, ik )
    !
    deltah = CMPLX(0.D0, 0.D0, kind=DP)
    DO ib = 1, nbnd
      DO jb = ib, nbnd
        !
        DO iwann = 1, num_wann
          IF ( gamma_only ) THEN 
            overlap_in = 2.D0 * SUM(DBLE(CONJG(evc(1:npw*npol,ib))*(evc0(1:npw*npol,iwann))))
            overlap_nj = 2.D0 * SUM(DBLE(CONJG(evc0(1:npw*npol,iwann))*(evc(1:npw*npol,jb))))
            IF (gstart == 2) THEN
               overlap_in = overlap_in - DBLE(CONJG(evc(1,ib))*(evc0(1,iwann)))
               overlap_nj = overlap_nj - DBLE(CONJG(evc0(1,iwann))*(evc(1,jb)))
            ENDIF
          ELSE
            overlap_in = SUM(CONJG(evc(1:npw*npol,ib))*(evc0(1:npw*npol,iwann)))
            overlap_nj = SUM(CONJG(evc0(1:npw*npol,iwann))*(evc(1:npw*npol,jb)))
          ENDIF
          CALL mp_sum (overlap_in, intra_bgrp_comm)
          CALL mp_sum (overlap_nj, intra_bgrp_comm)
          overlap = overlap_in*overlap_nj
          deltah(ib,jb) = deltah(ib,jb) + (0.5D0 - occ_mat(iwann)) * delta(iwann) * overlap
          !WRITE(*,'(3X, 2I5, 2F20.12, F20.12, 2F20.12)') ibnd, iwann, delta(iwann), occ_mat(iwann), overlap 
        ENDDO
        IF (ib /= jb) deltah(jb,ib) = CONJG(deltah(ib,jb))
        !
      ENDDO
    ENDDO
    !
    !
  END SUBROUTINE dki_hamiltonian
  !
  ! !----------------------------------------------------------------
  SUBROUTINE occupations (occ_mat)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kcw,           ONLY : num_wann, spin_component
    USE wvfct,                 ONLY : nbnd
    USE lsda_mod,              ONLY : nspin, lsda, isk
    USE klist,                 ONLY : ngk, nkstot, nks
    USE wvfct,                 ONLY : wg
    USE mp_bands,              ONLY : intra_bgrp_comm
    USE mp,                    ONLY : mp_sum
    USE mp_pools,              ONLY : inter_pool_comm
    USE control_flags,         ONLY : gamma_only
    USE gvect,                 ONLY : gstart
    !
    REAL(DP), INTENT(INOUT) :: occ_mat(num_wann)
    INTEGER :: iwann, ik, ibnd, ik_loc, spin_deg
    !
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    ! The canonical occupation matrix (fermi dirac or alike)
    occ_mat = REAL(0.D0, kind=DP)
    !
    ! ... P_n = \sum_kv f_kv <u_kv|w_kn><w_kn|u_kn> is a sum over ALL k: each pool
    ! accumulates the contribution of the k-points it owns (it can only read their
    ! KS orbitals from its own iuwfc buffer) and the partial sums are reduced over
    ! inter_pool_comm at the end.
    !
    DO iwann= 1, num_wann
      !
      DO ik_loc = 1, nks
        !
        IF ( lsda .AND. isk(ik_loc) /= spin_component ) CYCLE
        !
        ik = global_kpoint_index (nkstot, ik_loc) - (spin_component-1)*nkstot_eff
        !
        npw = ngk(ik_loc)
        lrwannfc = num_wann*npwx*npol
        CALL get_buffer ( evc0, lrwannfc, iuwfc_wann_allk, ik )
        CALL get_buffer ( evc, nwordwfc, iuwfc, ik_loc )
        !
        spin_deg=1
        IF(nspin == 1) spin_deg = 2
        overlap = CMPLX(0.D0, 0.D0, kind=DP)
        DO ibnd = 1, nbnd
          IF (gamma_only) THEN
             overlap = 2.D0*SUM(DBLE(CONJG(evc(1:npw*npol,ibnd))*(evc0(1:npw*npol,iwann))))
             IF (gstart ==2 ) overlap = overlap-1.D0*(DBLE(CONJG(evc(1,ibnd))*(evc0(1,iwann))))
          ELSE
             overlap = SUM(CONJG(evc(1:npw*npol,ibnd))*(evc0(1:npw*npol,iwann)))
          ENDIF
          CALL mp_sum (overlap, intra_bgrp_comm)
          overlap = CONJG(overlap)*overlap
          ! NB: wg is pool-local, so it takes the LOCAL index. It used to be indexed
          ! with the effective one, which also picked the wrong spin channel's weights
          ! for spin_component=2 (and coincided with the local index only for
          ! spin_component=1 without pools).
          occ_mat(iwann) = occ_mat(iwann) + wg(ibnd,ik_loc)/spin_deg * REAL(overlap)
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum ( occ_mat, inter_pool_comm )
    !
  END SUBROUTINE occupations

END SUBROUTINE koopmans_ham_proj
