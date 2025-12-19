!
! Copyright (C) 2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE diag_dense
  !----------------------------------------------------------------------------
  !!
  !! Direct diagonalization of the dense Hamiltonian matrix H(G,G').
  !! First build the matrix in the plane-wave basis, then diagonalize it using
  !! LAXlib routines.
  !!
  !! H(G,G') = δ_GG' × (k+G)² + V_loc(G-G') + Σ_ij |β_i(G)⟩ D_ij ⟨β_j(G')|
  !!
  !! RESTRICTIONS: NCPP only, no special features, serial only
  !!              (see diag_dense_check_compat)
  !!
  !! Inspired by the ParaBands code in BerkeleyGW
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: diag_dense_check_compat
  PUBLIC :: diag_dense_run_k
  PUBLIC :: diag_dense_test_hamiltonian
  !
  ! Derived type for Miller index lookup table
  !
  TYPE mill_lookup_type
     INTEGER, ALLOCATABLE :: map(:,:,:)
     !! 3D lookup array: map(nx,ny,nz) = G-vector index or 0 if not present
     INTEGER :: imin, imax
     !! Bounds in first dimension
     INTEGER :: jmin, jmax
     !! Bounds in second dimension
     INTEGER :: kmin, kmax
     !! Bounds in third dimension
  END TYPE mill_lookup_type
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  INTEGER FUNCTION lookup_mill_index( lookup, mill_vec ) RESULT(ig)
    !-----------------------------------------------------------------------
    !! Safe lookup of G-vector index from Miller indices
    !! Returns 0 if Miller indices are out of bounds or G-vector not present
    !
    IMPLICIT NONE
    !
    TYPE(mill_lookup_type), INTENT(IN) :: lookup
    !! Lookup table structure
    INTEGER, INTENT(IN) :: mill_vec(3)
    !! Miller indices to look up
    !
    ! Check bounds
    IF ( mill_vec(1) < lookup%imin .OR. mill_vec(1) > lookup%imax .OR. &
         mill_vec(2) < lookup%jmin .OR. mill_vec(2) > lookup%jmax .OR. &
         mill_vec(3) < lookup%kmin .OR. mill_vec(3) > lookup%kmax ) THEN
       ig = 0
    ELSE
       ! Lookup in table
       ig = lookup%map(mill_vec(1), mill_vec(2), mill_vec(3))
    ENDIF
    !
  END FUNCTION lookup_mill_index
  !
  !-----------------------------------------------------------------------
  SUBROUTINE build_mill_lookup_table( lookup )
    !-----------------------------------------------------------------------
    !! Build 3D lookup table: Miller indices (Gx,Gy,Gz) → global G-vector index
    !!
    !! Uses GLOBAL Miller indices to ensure all G-vectors are indexed,
    !! including those from other processors in parallel execution.
    !!
    !! Returns global G-vector indices (ig_l2g) for use across processors.
    !
    USE gvect,         ONLY : mill, ngm, ig_l2g
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_max, mp_min, mp_sum
    !
    IMPLICIT NONE
    !
    TYPE(mill_lookup_type), INTENT(OUT) :: lookup
    !! Lookup table structure (output)
    !
    INTEGER :: ig, ierr
    INTEGER :: mill_max(3), mill_min(3)
    REAL(DP) :: mem_mb
    !
    ! Determine bounds from actual Miller indices using intrinsic functions
    !
    mill_max(1) = MAXVAL( mill(1, 1:ngm) )
    mill_max(2) = MAXVAL( mill(2, 1:ngm) )
    mill_max(3) = MAXVAL( mill(3, 1:ngm) )
    !
    mill_min(1) = MINVAL( mill(1, 1:ngm) )
    mill_min(2) = MINVAL( mill(2, 1:ngm) )
    mill_min(3) = MINVAL( mill(3, 1:ngm) )
    !
    ! Get global maximum/minimum across all processors
    CALL mp_max( mill_max, intra_bgrp_comm )
    CALL mp_min( mill_min, intra_bgrp_comm )
    !
    ! Set bounds (no padding needed)
    lookup%imin = mill_min(1)
    lookup%imax = mill_max(1)
    lookup%jmin = mill_min(2)
    lookup%jmax = mill_max(2)
    lookup%kmin = mill_min(3)
    lookup%kmax = mill_max(3)
    !
    ! Allocate 3D array
    ALLOCATE(lookup%map(lookup%imin:lookup%imax, &
                        lookup%jmin:lookup%jmax, &
                        lookup%kmin:lookup%kmax))
    lookup%map(:,:,:) = 0
    !
    ! Fill lookup table for the global G-vector indices.
    ! Each processor fills its own G-vectors, and then we sum across all.
    !
    DO ig = 1, ngm
       lookup%map(mill(1,ig), mill(2,ig), mill(3,ig)) = ig_l2g(ig)
    ENDDO
    !
    CALL mp_sum(lookup%map, intra_bgrp_comm)
    !
  END SUBROUTINE build_mill_lookup_table
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_check_compat()
    !-----------------------------------------------------------------------
    !! Check if calculation settings are compatible with full-band diagonalization
    !!
    !! Calls errore and stops if incompatible
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE control_flags,    ONLY : gamma_only, use_para_diag
    USE realus,           ONLY : real_space
    USE noncollin_module, ONLY : noncolin
    USE uspp,             ONLY : okvan
    USE ldaU,             ONLY : lda_plus_u
    USE bp,               ONLY : lelfield
    USE xc_lib,           ONLY : exx_is_active, xclib_dft_is
    USE fft_base,         ONLY : dffts
    USE wvfct,            ONLY : npwx
    !
    IMPLICIT NONE
    !
    REAL(DP) :: mem_gb
    !
    ! Check for incompatible features
    !
    IF ( okvan ) CALL errore('diag_dense_check_compat', &
       'Full-band diagonalization requires NCPP (no USPP/PAW)', 1)
    !
    IF ( gamma_only ) CALL errore('diag_dense_check_compat', &
       'gamma_only not supported in full-band diagonalization', 1)
    !
    IF ( real_space ) CALL errore('diag_dense_check_compat', &
       'real_space algorithms not supported in full-band diagonalization', 1)
    !
    IF ( noncolin ) CALL errore('diag_dense_check_compat', &
       'non-collinear magnetism not supported in full-band diagonalization', 1)
    !
    IF ( xclib_dft_is('meta') ) CALL errore('diag_dense_check_compat', &
       'meta-GGA functionals not supported in full-band diagonalization', 1)
    !
    IF ( lda_plus_u ) CALL errore('diag_dense_check_compat', &
       'DFT+U not supported in full-band diagonalization', 1)
    !
    IF ( exx_is_active() ) CALL errore('diag_dense_check_compat', &
       'Exact exchange/hybrid functionals not supported in full-band diagonalization', 1)
    !
    IF ( lelfield ) CALL errore('diag_dense_check_compat', &
       'Electric field (Berry phase) not supported in full-band diagonalization', 1)
    !
    IF ( dffts%has_task_groups ) CALL errore('diag_dense_check_compat', &
       'Task groups not supported in full-band diagonalization', 1)
    !
    ! Memory estimate for full-band diagonalization
    ! FIXME: Call it somewhere else when npwx is initialized
    !
    mem_gb = REAL(npwx,DP)**2 * 16.0_DP / 1024.0_DP**3
    WRITE(stdout,'(/,5X,"Full-band diagonalization memory estimate:")')
    WRITE(stdout,'(5X,"npwx =",I8,", H matrix requires",F8.2," GB")') npwx, mem_gb
    !
    IF ( mem_gb > 1.0_DP ) THEN
       WRITE(stdout,'(5X,"WARNING: Large memory requirement (>1 GB)")')
       WRITE(stdout,'(5X,"         Diagonalization may take several minutes")')
    ENDIF
    !
  END SUBROUTINE diag_dense_check_compat
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_hamiltonian_kinetic( hmat, ik, npw )
    !-----------------------------------------------------------------------
    !! Build kinetic energy contribution (diagonal)
    !
    USE kinds,  ONLY : DP
    USE wvfct,  ONLY : g2kin
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(npw,npw)
    !! Hamiltonian matrix
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    !
    INTEGER :: ig
    !
    CALL start_clock( 'ham_kin' )
    !
    DO ig = 1, npw
       hmat(ig,ig) = CMPLX(g2kin(ig), 0.0_DP, KIND=DP)
    ENDDO
    !
    CALL stop_clock( 'ham_kin' )
    !
  END SUBROUTINE diag_dense_hamiltonian_kinetic
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_hamiltonian_vloc( hmat, ik, npw )
    !-----------------------------------------------------------------------
    !! Build local potential contribution V_loc(G-G')
    !
    USE kinds,         ONLY : DP
    USE klist,         ONLY : igk_k
    USE gvect,         ONLY : mill, ngm, ngm_g, ig_l2g
    USE lsda_mod,      ONLY : current_spin
    USE scf,           ONLY : vrs
    USE fft_base,      ONLY : dffts
    USE fft_interfaces,ONLY : fwfft
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(npw,npw)
    !! Hamiltonian matrix
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    !
    ! Local variables
    COMPLEX(DP), ALLOCATABLE :: vrs_g(:)
    !! Local potential in G-space
    COMPLEX(DP), ALLOCATABLE :: vrs_aux(:)
    !! FFT workspace for V_loc
    TYPE(mill_lookup_type) :: mill_lookup
    !! Lookup table structure for Miller indices → G-vector index
    INTEGER :: ig, jg, igg
    INTEGER :: mill_diff(3)
    !
    CALL start_clock( 'ham_vloc' )
    !
    ALLOCATE(vrs_g(ngm_g))
    ALLOCATE(vrs_aux(dffts%nnr))
    vrs_g(:) = (0.0_DP, 0.0_DP)
    !
    ! Transform V_loc to G-space
    vrs_aux(:) = CMPLX(vrs(:,current_spin), 0.0_DP, KIND=DP)
    CALL fwfft('Rho', vrs_aux, dffts)
    !
    ! Extract G-space components using global G-vector indices
    DO ig = 1, ngm
       vrs_g(ig_l2g(ig)) = vrs_aux(dffts%nl(ig))
    ENDDO
    !
    ! Sum across all processors
    CALL mp_sum( vrs_g, intra_bgrp_comm )
    !
    ! Build Miller index lookup table
    CALL build_mill_lookup_table( mill_lookup )
    !
    ! Build V_loc contributions (off-diagonal)
    DO jg = 1, npw
       DO ig = 1, npw
          ! Compute Miller indices difference G - G'
          mill_diff(1) = mill(1,igk_k(ig,ik)) - mill(1,igk_k(jg,ik))
          mill_diff(2) = mill(2,igk_k(ig,ik)) - mill(2,igk_k(jg,ik))
          mill_diff(3) = mill(3,igk_k(ig,ik)) - mill(3,igk_k(jg,ik))
          !
          ! O(1) lookup using safe function
          igg = lookup_mill_index( mill_lookup, mill_diff )
          IF ( igg > 0 ) THEN
             hmat(ig,jg) = hmat(ig,jg) + vrs_g(igg)
          ENDIF
       ENDDO
    ENDDO
    !
    ! Cleanup
    DEALLOCATE( vrs_g, vrs_aux, mill_lookup%map )
    !
    CALL stop_clock( 'ham_vloc' )
    !
  END SUBROUTINE diag_dense_hamiltonian_vloc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_hamiltonian_vnl( hmat, ik, npw )
    !-----------------------------------------------------------------------
    !! Build nonlocal pseudopotential contribution
    !! H(G,G') += Σ_ij |β_i(G)⟩ D_ij ⟨β_j(G')|
    !
    USE kinds,      ONLY : DP
    USE wvfct,      ONLY : npwx
    USE lsda_mod,   ONLY : current_spin
    USE uspp,       ONLY : nkb, vkb, deeq, ofsbeta
    USE uspp_param, ONLY : nh, nhm
    USE ions_base,  ONLY : nat, ntyp => nsp, ityp
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(npw,npw)
    !! Hamiltonian matrix
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    !
    ! Local variables
    COMPLEX(DP), ALLOCATABLE :: ps(:,:)
    !! Temp array for KB projector application
    COMPLEX(DP), ALLOCATABLE :: deeaux(:,:)
    !! Complex copy of deeq for ZGEMM
    INTEGER :: nt, na, ib, jb, jg
    !
    CALL start_clock( 'ham_vnl' )
    !
    IF ( nkb > 0 ) THEN
       !
       ALLOCATE(ps(nkb, npw))
       ps(:,:) = (0.0_DP, 0.0_DP)
       !
       ! Loop over atom types and atoms
       DO nt = 1, ntyp
          IF ( nh(nt) == 0 ) CYCLE
          !
          ALLOCATE(deeaux(nh(nt), nh(nt)))
          !
          DO na = 1, nat
             IF ( ityp(na) /= nt ) CYCLE
             !
             ! Copy real deeq into complex deeaux
             deeaux(1:nh(nt), 1:nh(nt)) = CMPLX(deeq(1:nh(nt), 1:nh(nt), na, current_spin), &
                                                 0.0_DP, KIND=DP)
             !
             ! For each G' (column jg): compute ps(beta) = D * <beta|G'>
             DO jg = 1, npw
                DO ib = 1, nh(nt)
                   ps(ofsbeta(na)+ib, jg) = (0.0_DP, 0.0_DP)
                   DO jb = 1, nh(nt)
                      ps(ofsbeta(na)+ib, jg) = ps(ofsbeta(na)+ib, jg) + &
                           deeaux(ib,jb) * CONJG(vkb(jg, ofsbeta(na)+jb))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO  ! na
          !
          DEALLOCATE( deeaux )
       ENDDO  ! nt
       !
       ! Add the nonlocal contribution: H(ig,jg) += Σ_beta vkb(ig,beta) * ps(beta,jg)
       CALL ZGEMM( 'N', 'N', npw, npw, nkb, (1.0_DP, 0.0_DP), &
                   vkb, npwx, ps, nkb, (1.0_DP, 0.0_DP), hmat, npw )
       !
       DEALLOCATE( ps )
       !
    ENDIF
    !
    CALL stop_clock( 'ham_vnl' )
    !
  END SUBROUTINE diag_dense_hamiltonian_vnl
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_check_herm( hmat, npw )
    !-----------------------------------------------------------------------
    !! Check and report Hermiticity of the Hamiltonian matrix
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: hmat(npw,npw)
    !! Hamiltonian matrix
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    !
    REAL(DP) :: max_herm_err
    INTEGER :: ig, jg
    !
    max_herm_err = 0.0_DP
    DO jg = 1, MIN(npw, 100)  ! Check first 100 for speed
       DO ig = jg+1, MIN(npw, 100)
          max_herm_err = MAX(max_herm_err, ABS(hmat(ig,jg) - CONJG(hmat(jg,ig))))
       ENDDO
    ENDDO
    !
    IF ( max_herm_err > 1.0e-6_DP ) THEN
       WRITE(stdout,'(5X,"WARNING: H matrix not Hermitian, max error =",E12.4)') max_herm_err
    ELSE
       WRITE(stdout,'(5X,"Hermiticity check: max error =",E12.4," [OK]")') max_herm_err
    ENDIF
    !
  END SUBROUTINE diag_dense_check_herm
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_diag( hmat, smat, npw, nbnd, evc, et_k, notconv )
    !-----------------------------------------------------------------------
    !! Diagonalize the Hamiltonian matrix using LAXlib
    !
    USE kinds,     ONLY : DP
    USE wvfct,     ONLY : npwx
    USE io_global, ONLY : stdout
    USE mp_bands,  ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
    USE LAXlib,    ONLY : diaghg
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: hmat(npw,npw)
    !! Hamiltonian matrix
    COMPLEX(DP), INTENT(INOUT) :: smat(npw,npw)
    !! Overlap matrix
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    INTEGER, INTENT(IN) :: nbnd
    !! number of bands
    COMPLEX(DP), INTENT(OUT) :: evc(npwx,nbnd)
    !! eigenvectors
    REAL(DP), INTENT(OUT) :: et_k(nbnd)
    !! eigenvalues
    INTEGER, INTENT(OUT) :: notconv
    !! number of non-converged bands (always 0 for direct method)
    !
    REAL(DP), ALLOCATABLE :: et_tmp(:)
    !! Temporary array for eigenvalues
    COMPLEX(DP), ALLOCATABLE :: evc_tmp(:, :)
    !! Temporary array for eigenvectors
    !
    ALLOCATE(et_tmp(npw))
    ALLOCATE(evc_tmp(npw, nbnd))
    !
    WRITE(stdout,'(5X,"Diagonalizing",I5,"x",I5," matrix...")') npw, npw
    !
    CALL start_clock( 'ham_diag' )
    CALL diaghg( npw, nbnd, hmat, smat, npw, et_tmp, evc_tmp, &
                 me_bgrp, root_bgrp, intra_bgrp_comm )
    CALL stop_clock( 'ham_diag' )
    !
    evc(:, :) = (0.0_DP, 0.0_DP)
    et_k(:) = et_tmp(1:nbnd)
    evc(1:npw, 1:nbnd) = evc_tmp(1:npw, 1:nbnd)
    !
    notconv = 0  ! Always converged (direct diagonalization)
    !
  END SUBROUTINE diag_dense_diag
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_run_k( ik, npw, nbnd, evc, et_k, notconv )
    !-----------------------------------------------------------------------
    !! Construct and diagonalize explicit Hamiltonian matrix for k-point ik
    !!
    !! This is the main coordinator that:
    !! 1. Allocates H and S matrices
    !! 2. Builds identity S matrix (for NCPP)
    !! 3. Calls hamiltonian_kinetic, hamiltonian_vloc, hamiltonian_vnl
    !! 4. Checks hermiticity
    !! 5. Calls diagonalize_hamiltonian
    !! 6. Verifies against h_psi
    !
    USE kinds,     ONLY : DP
    USE wvfct,     ONLY : npwx
    USE io_global, ONLY : stdout
    USE uspp,      ONLY : nkb
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    INTEGER, INTENT(IN) :: nbnd
    !! number of bands
    COMPLEX(DP), INTENT(OUT) :: evc(npwx,nbnd)
    !! eigenvectors (output)
    REAL(DP), INTENT(OUT) :: et_k(nbnd)
    !! eigenvalues (output)
    INTEGER, INTENT(OUT) :: notconv
    !! number of non-converged bands (always 0 for direct method)
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: hmat(:,:)
    !! Hamiltonian matrix (npw, npw)
    COMPLEX(DP), ALLOCATABLE :: smat(:,:)
    !! Overlap matrix (identity for NCPP)
    INTEGER :: ierr, ig
    !! allocation error code, loop index
    REAL(DP) :: mem_gb
    !! memory estimate
    !
    CALL start_clock( 'hamiltonian_dense' )
    !
    ! Report memory requirements
    !
    mem_gb = REAL(npw,DP)**2 * 16.0_DP / 1024.0_DP**3
    WRITE(stdout,'(/,5X,"Dense H construction for k-point",I5)') ik
    WRITE(stdout,'(5X,"npw =",I8,", Matrix memory:",F8.3," GB")') npw, mem_gb
    !
    ! Allocate matrices
    !
    ALLOCATE( hmat(npw, npw), STAT=ierr )
    IF ( ierr /= 0 ) THEN
       WRITE(stdout,'(5X,"ERROR: Cannot allocate H matrix, need",F8.2," GB")') mem_gb
       CALL errore('construct_hamiltonian_k', 'H matrix allocation failed', ierr)
    ENDIF
    !
    ALLOCATE( smat(npw, npw), STAT=ierr )
    IF ( ierr /= 0 ) CALL errore('construct_hamiltonian_k', 'S matrix allocation failed', ierr)
    !
    ! Initialize matrices
    !
    hmat(:,:) = (0.0_DP, 0.0_DP)
    smat(:,:) = (0.0_DP, 0.0_DP)
    !
    ! Build identity overlap matrix (for NCPP)
    !
    DO ig = 1, npw
       smat(ig,ig) = (1.0_DP, 0.0_DP)
    ENDDO
    !
    ! Build Hamiltonian terms
    !
    CALL diag_dense_hamiltonian_kinetic(hmat, ik, npw)
    CALL diag_dense_hamiltonian_vloc(hmat, ik, npw)
    CALL diag_dense_hamiltonian_vnl(hmat, ik, npw)
    !
    ! Check Hermiticity
    !
    CALL diag_dense_check_herm(hmat, npw)
    !
    ! Verify H matrix against h_psi (always run verification for testing)
    !
    CALL diag_dense_test_hamiltonian(hmat, ik, npw, 1.0e-8_DP)
    !
    ! Diagonalize
    !
    CALL diag_dense_diag(hmat, smat, npw, nbnd, evc, et_k, notconv)
    !
    ! Deallocate
    !
    DEALLOCATE(hmat)
    DEALLOCATE(smat)
    !
    CALL stop_clock( 'hamiltonian_dense' )
    !
    ! Print timing breakdown
    !
    WRITE(stdout,'(/,5X,"Timing breakdown for dense H construction:")')
    CALL print_clock( 'ham_kin' )
    CALL print_clock( 'ham_vloc' )
    CALL print_clock( 'ham_vnl' )
    CALL print_clock( 'ham_diag' )
    CALL print_clock( 'hamiltonian_dense' )
    !
  END SUBROUTINE diag_dense_run_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_dense_test_hamiltonian( hmat, ik, npw, tolerance )
    !-----------------------------------------------------------------------
    !! Verify constructed H matrix by comparing H*v with h_psi(v)
    !! for a random test vector v
    !
    USE kinds,    ONLY : DP
    USE wvfct,    ONLY : npwx
    USE noncollin_module, ONLY : npol
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    COMPLEX(DP), INTENT(IN) :: hmat(npw,npw)
    !! Hamiltonian matrix to verify
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    INTEGER, INTENT(IN) :: npw
    !! number of plane waves
    REAL(DP), INTENT(IN) :: tolerance
    !! tolerance for passing test
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: v(:)
    !! random test vector
    COMPLEX(DP), ALLOCATABLE :: hv_mat(:)
    !! H*v from matrix multiply
    COMPLEX(DP), ALLOCATABLE :: hv_hpsi(:)
    !! H*v from h_psi
    REAL(DP), ALLOCATABLE :: rand_vec(:)
    !! random numbers
    REAL(DP) :: vnorm, max_error, rms_error
    INTEGER :: ig, ierr
    !
    EXTERNAL :: h_psi
    ! subroutine h_psi(lda,n,m,psi,hpsi)
    !
    WRITE(stdout,'(/,5X,"Verifying H matrix against h_psi...")')
    !
    ! Allocate arrays
    !
    ALLOCATE( v(npw), hv_mat(npw), hv_hpsi(npwx), rand_vec(npw), STAT=ierr )
    IF ( ierr /= 0 ) CALL errore('verify_ham_vs_hpsi', 'allocation failed', ierr)
    !
    ! Generate random normalized vector
    !
    CALL RANDOM_NUMBER(rand_vec)
    v(:) = CMPLX(rand_vec(:), 0.0_DP, KIND=DP)
    vnorm = SQRT(SUM(ABS(v)**2))
    v(:) = v(:) / vnorm
    !
    ! Compute H*v using matrix multiply
    !
    hv_mat(:) = (0.0_DP, 0.0_DP)
    DO ig = 1, npw
       hv_mat(ig) = SUM(hmat(ig,:) * v(:))
    ENDDO
    !
    ! Compute H*v using h_psi
    !
    hv_hpsi(:) = (0.0_DP, 0.0_DP)
    CALL h_psi( npwx, npw, 1, v, hv_hpsi )
    !
    ! Compare results
    !
    max_error = 0.0_DP
    rms_error = 0.0_DP
    DO ig = 1, npw
       max_error = MAX(max_error, ABS(hv_mat(ig) - hv_hpsi(ig)))
       rms_error = rms_error + ABS(hv_mat(ig) - hv_hpsi(ig))**2
    ENDDO
    rms_error = SQRT(rms_error / REAL(npw,DP))
    !
    ! Report results
    !
    WRITE(stdout,'(5X,"max |H_mat*v - h_psi(v)| =",E12.4)') max_error
    WRITE(stdout,'(5X,"RMS |H_mat*v - h_psi(v)| =",E12.4)') rms_error
    !
    IF ( max_error < tolerance ) THEN
       WRITE(stdout,'(5X,"Verification: PASSED")')
    ELSE
       WRITE(stdout,'(5X,"WARNING: Verification FAILED!")')
       WRITE(stdout,'(5X,"         Expected max error <",E12.4)') tolerance
    ENDIF
    !
    DEALLOCATE( v, hv_mat, hv_hpsi, rand_vec )
    !
  END SUBROUTINE diag_dense_test_hamiltonian
  !
END MODULE diag_dense
