!
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
SUBROUTINE koopmans_ham_uniq ( dH_wann )
  !---------------------------------------------------------------------
  !
  ! Here the KI hamiltonian is written in terms of projectors on Wannnier
  ! functions:
  ! \Delta H_KI = \sum_nm |w_n> \Delta H_nm < w_m| where \Delta H_nm = <w_n | h_m |w_m> 
  ! computed in the standard way (see dH_ki_wann.f90). 
  ! Then we build and Diagonalize the KI hamiltoinan H_KS+\Delta H_KI on the basis of the 
  ! KS orbitals from the NSCF calculation: espilon_i = Diag [ < \phi_i | H_KS + \Delta H_KI | phi_j > ] 
  ! < \phi_i | H_KS | phi_j > = \delta_ij \epsilon_i^KS 
  ! < \phi_i | \Delta H_KI | phi_j > = \sum_nm <phi_i|w_n> \Delta H_nm <w_m|phi_j> 
  !
  ! NB: In principle one can use any other basis or iterative digonalization technique.  
  
  USE io_global,             ONLY : stdout, ionode
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : xk, ngk, nkstot, nks
  USE lsda_mod,              ONLY : lsda, isk
  USE control_kcw,           ONLY : num_wann, evc0, spin_component, &
                                    num_wann_occ, iuwfc_wann_allk, nkstot_eff, &
                                    kcw_iverbosity
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx, npw, et, nbnd
  USE units_lr,              ONLY : iuwfc
  USE wavefunctions,         ONLY : evc
  USE buffers,               ONLY : get_buffer, save_buffer
  !
  USE io_files,              ONLY : nwordwfc
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum, mp_max, mp_min
  USE mp_pools,              ONLY : inter_pool_comm
  USE noncollin_module,      ONLY : npol
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !! The global index of a local (pool) k-point
  !
  ! ik_loc is the LOCAL (pool) k-point index: the KS orbitals (iuwfc), the number
  ! of PWs (ngk) and the eigenvalues (et) are all pool-local arrays, so they must
  ! be addressed with it. ik is the "effective" (1:nkstot_eff) index of the current
  ! spin channel, used for the pool-replicated quantities: dH_wann and the ALL-k
  ! Wannier-gauge buffer iuwfc_wann_allk.
  INTEGER :: ik, ik_loc
  !
  ! Per-k results, stashed at the effective k index and gathered across pools after
  ! the loop so that the full, correctly ordered table is printed once from ionode.
  ! Printing from inside the loop would only ever reach the log for the k-points
  ! owned by ionode's own pool.
  REAL(DP), ALLOCATABLE :: eigvl_ks_all(:,:), eigvl_ki_all(:,:), eigvl_pert_all(:,:)
  REAL(DP), ALLOCATABLE :: xk_all(:,:)
  ! xk is pool-local too, so the k coordinates used as print labels are gathered as well
  REAL(DP), ALLOCATABLE :: ki_spec(:,:,:)
  ! kcw_iverbosity>1 only: the empty-state spectrum as a function of the size of the
  ! diagonalized subspace; at most the first 10 eigenvalues of each are printed
  !
  ! the KI hamiltonian on the Wannier basis <w_i|dh_j|w_j> 
  COMPLEX(DP), INTENT (IN) :: dH_wann(nkstot_eff,num_wann,num_wann)
  COMPLEX(DP), ALLOCATABLE :: dH_wann_aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_aux(:,:)
  ! 
  ! the KI operator on the KS basis of the NSCF calculation
  COMPLEX(DP) :: deltah(nbnd,nbnd)
  ! the new hamitonain, and the new eigenvalues and eigenvectors at a given k-point
  COMPLEX(DP) :: ham(nbnd,nbnd), eigvc(nbnd,nbnd)
  !
  ! The new eigenalues 
  REAL(DP) :: eigvl(nbnd)
  REAL(DP) :: eigvl_pert(nbnd)
  REAL(DP) :: eigvl_ks(nbnd)
  !
  INTEGER :: i, ibnd
  ! 
  REAL(DP) :: ehomo, elumo
  REAL(DP) :: ehomo_ks, elumo_ks
  REAL(DP) :: ehomo_pert, elumo_pert
  REAL(DP), EXTERNAL :: get_clock
  EXTERNAL :: ZGEMM, CDIAGH
  !
  COMPLEX(DP), ALLOCATABLE :: ham_aux(:,:)
  REAL(DP), ALLOCATABLE :: eigvl_ki(:)
  COMPLEX(DP), ALLOCATABLE :: eigvc_ki(:,:)
  INTEGER i_start, i_end, k
  !
  !
  ehomo=-1D+6
  elumo=+1D+6
  ehomo_ks=-1D+6
  elumo_ks=+1D+6
  ehomo_pert=-1D+6
  elumo_pert=+1D+6
  !
  WRITE( stdout, '(/,5X, "INFO: BUILD and DIAGONALIZE the KI HAMILTONIAN")')
  WRITE( stdout, '(  5X, "INFO: Unique Hamiltonian scheme using projectors:")')
  WRITE( stdout, '(  5X, "      build and diagonalize the KI Hamiltonian in")')
  WRITE( stdout, '(  5X, "      the basis of KS orbitals")')
  !
  ALLOCATE ( dH_wann_aux(num_wann, num_wann) )
  ALLOCATE ( evc_aux(npwx*npol, nbnd) )
  ALLOCATE ( eigvl_ks_all(nbnd, nkstot_eff), eigvl_ki_all(nbnd, nkstot_eff), &
             eigvl_pert_all(nbnd, nkstot_eff), xk_all(3, nkstot_eff) )
  eigvl_ks_all = 0.D0; eigvl_ki_all = 0.D0; eigvl_pert_all = 0.D0; xk_all = 0.D0
  IF (kcw_iverbosity .gt. 1 .AND. nbnd > num_wann_occ) THEN
    ALLOCATE ( ki_spec(10, nbnd-num_wann_occ, nkstot_eff) )
    ki_spec = 0.D0
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
    dH_wann_aux(:,:) = dH_wann(ik, :,:)
    !
    ! Unique Hamiltonian diagonalized on the KS basis of the NSF calculation
    !
    CALL get_buffer ( evc, nwordwfc, iuwfc, ik_loc )
    npw = ngk(ik_loc)
    xk_all(:,ik) = xk(:,ik_loc)
    !
    ! The KS Hamiltonian in the KS basis
    ham(:,:)=CMPLX(0.D0, 0.D0, kind=DP)
    DO i = 1, nbnd
      ham(i,i)    = et(i,ik_loc)
      eigvl_ks(i) = et(i,ik_loc)
    ENDDO
    !
    ehomo_ks = MAX ( ehomo_ks, eigvl_ks(num_wann_occ  ) )
    IF (nbnd > num_wann_occ) elumo_ks = MIN ( elumo_ks, eigvl_ks(num_wann_occ+1) )
    !
    ! The Delta H_KI_ij = \sum_nm <phi_i|w_n> \Delta H_nm <w_m|phi_j>
    CALL dki_hamiltonian (evc, ik, nbnd, dH_wann_aux(:,:), deltah)
    !
    ! Add to the KS Hamiltonian
    ham(:,:) = ham(:,:) + deltah(:,:)
    !
#ifdef DEBUG
    WRITE(stdout, '(/, "dKI Hamiltonian at k = ", i4)') ik
    DO k = 1, 10
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(deltah(k,i)),AIMAG(deltah(k,i)), i=1,10)
    ENDDO
    !
    WRITE(stdout, '(/, "KI Hamiltonian at k = ", i4)') ik
    DO k = 1, 10
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(ham(k,i)),AIMAG(ham(k,i)), i=1,10)
    ENDDO
#endif
    !
    ! Because we have defined a uniq KI Hamiltonian, we can do a perturbative approach
    ! i.e. we keep only the diagonal part of the KI Hamiltoniana
    DO i = 1, nbnd
      eigvl_pert(i) = et(i,ik_loc) + DBLE(deltah(i,i))
    ENDDO
    ehomo_pert = MAX ( ehomo_pert, eigvl_pert(num_wann_occ ) )
    IF (nbnd > num_wann_occ) elumo_pert = MIN ( elumo_pert, eigvl_pert(num_wann_occ+1 ) )
    !
    IF (kcw_iverbosity .gt. 1 .AND. nbnd > num_wann_occ ) THEN
      !
      ! Stash the empty-state spectrum; printed after the loop (see below)
      DO k = 1, nbnd-num_wann_occ
         !
         i_start = num_wann_occ+1; i_end = num_wann_occ+k
         !
         ALLOCATE (ham_aux(k,k), eigvl_ki(k), eigvc_ki(k,k))
         ham_aux(1:k,1:k) = ham(i_start:i_end,i_start:i_end)
         !
         CALL cdiagh( k, ham_aux, k, eigvl_ki, eigvc_ki )
         !
         ki_spec(1:MIN(k,10), k, ik) = eigvl_ki(1:MIN(k,10))   ! First 10 eigenvalues
         !
         DEALLOCATE (ham_aux)
         DEALLOCATE (eigvl_ki, eigvc_ki)
         !
      ENDDO
    ENDIF
    !
    ! Diagonalize the KI Hamiltonian
    CALL CDIAGH( nbnd, ham, nbnd, eigvl, eigvc )
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
  IF ( ALLOCATED(ki_spec) ) CALL mp_sum ( ki_spec, inter_pool_comm )
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
      IF ( ALLOCATED(ki_spec) ) THEN
        WRITE(stdout,'(8x, "INFO: Empty states spectrum as a function of the # of orbitals")')
        DO k = 1, nbnd-num_wann_occ
          WRITE(stdout,'(8x, I3, 10F10.4)') k, ki_spec(1:MIN(k,10), k, ik)*rytoev
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
  DEALLOCATE (dH_wann_aux)
  DEALLOCATE (evc_aux)
  DEALLOCATE (eigvl_ks_all, eigvl_ki_all, eigvl_pert_all, xk_all)
  IF ( ALLOCATED(ki_spec) ) DEALLOCATE (ki_spec)
  !
  9043 FORMAT(/,8x, 'KS  highest occupied level (ev): ',F10.4 )
  9042 FORMAT(/,8x, 'KS  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
  9045 FORMAT(  8x, 'KI  highest occupied level (ev): ',F10.4 )
  9044 FORMAT(  8x, 'KI  highest occupied, lowest unoccupied level (ev): ',2F10.4 )
  9047 FORMAT(  8x, 'pKI highest occupied level (ev): ',F10.4 )
  9046 FORMAT(  8x, 'pKI highest occupied, lowest unoccupied level (ev): ',2F10.4 )
  9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
  901 FORMAT('          total cpu time spent up to now is ',F10.1,' secs' )
  !
  RETURN
  CONTAINS
  !
  ! !----------------------------------------------------------------
  SUBROUTINE dki_hamiltonian (evc, ik, h_dim, delta, deltah)
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
    COMPLEX(DP) :: delta (num_wann,num_wann)
    COMPLEX(DP), INTENT(IN) :: evc(npwx*npol,h_dim)
    INTEGER, INTENT(IN) :: ik
    COMPLEX(DP), INTENT(OUT) :: deltah(h_dim,h_dim)
    !
    INTEGER :: lrwannfc, ib, jb, nwann, mwann
    COMPLEX(DP) :: overlap
    !
    COMPLEX (DP) :: overlap_mat(nbnd,num_wann)
    REAL(DP), ALLOCATABLE :: overlap_mat_real(:,:)
    !
    EXTERNAL :: ZGEMM
    !
    lrwannfc = num_wann*npwx*npol
    CALL get_buffer ( evc0, lrwannfc, iuwfc_wann_allk, ik )
    !
    deltah = CMPLX(0.D0, 0.D0, kind=DP)
    !
    IF (gamma_only) THEN 
      ALLOCATE ( overlap_mat_real(nbnd,num_wann) )
      CALL DGEMM( 'C', 'N', nbnd, num_wann, 2*npw, 2.0_DP, evc, 2*npwx, evc0, &
             2*npwx, 0.0_DP, overlap_mat_real, nbnd )
     IF ( gstart == 2 ) &
        CALL DGER( nbnd, num_wann, -1.0_DP, evc , 2*npwx, evc0, 2*npwx, overlap_mat_real, nbnd )
    ELSE
      CALL ZGEMM( 'C','N', nbnd, num_wann, npwx*npol, ONE, evc, npwx*npol, evc0, npwx*npol, &
      ZERO, overlap_mat, nbnd) 
    ENDIF
    !
    IF (gamma_only) THEN 
       overlap_mat = CMPLX(overlap_mat_real, 0.D0, kind=DP)
       DEALLOCATE ( overlap_mat_real )
    ENDIF
    !
    CALL mp_sum (overlap_mat, intra_bgrp_comm)
    !
    DO ib = 1, nbnd
      DO jb = ib, nbnd
        !
        DO nwann = 1, num_wann
          DO mwann = 1, num_wann
            ! 
            overlap = (overlap_mat(ib,nwann )) * CONJG(overlap_mat(jb,mwann))
            deltah(ib,jb) = deltah(ib,jb) + delta(nwann,mwann) * overlap
            !WRITE(*,'(5X, 2I5, 2F20.12, 2F20.12, 2F20.12)') nwann, mwann, delta(nwann,mwann), overlap, deltah(ib,jb)
          ENDDO
        ENDDO
        !WRITE(*,'(3X, 2I5, 2F20.12)') ib, jb, deltah(ib,jb)
        IF (ib /= jb) deltah(jb,ib) = CONJG(deltah(ib,jb))
        !
      ENDDO
    ENDDO
    !
    !
  END SUBROUTINE dki_hamiltonian
  !
  !
END SUBROUTINE koopmans_ham_uniq
