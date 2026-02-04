! Copyrigh(C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE exx2
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only, use_gpu
  USE exx_base,             ONLY : dfftt , exxbuff, exxbuff_d, npwt, x_nbnd_occ, &
                                   ibnd_start, ibnd_end
  USE noncollin_module,     ONLY : noncolin, npol
  !
  INTEGER :: ibnd_buff_start
  !! starting buffer index used in bgrp parallelization
  INTEGER :: ibnd_buff_end
  !! ending buffer index used in bgrp parallelization
  !
 CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit2( )
    !------------------------------------------------------------------------
    !! This subroutine is run before the first H_psi() of each iteration. 
    !! It saves the wavefunctions for the right density matrix, in real space.
    !
    USE wavefunctions,        ONLY : psic
    USE io_files,             ONLY : iunwfc_exx
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx
    USE klist,                ONLY : ngk, nks, nkstot
    USE mp_exx,               ONLY : me_egrp, init_index_over_band,  &
                                     inter_egrp_comm,           &
                                     iexx_start, iexx_end, &
                                     all_start, all_end
    USE mp,                   ONLY : mp_bcast
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : okvan
    USE us_exx,               ONLY : rotate_becxx
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_fock_kernel
    USE mp_orthopools,        ONLY : intra_orthopool_comm
    USE exx_base,             ONLY : nkqs, xkq_collect, index_xk, index_sym,  &
                                     rir, working_pool
    USE exx_band,             ONLY : change_data_structure, nwordwfc_exx, &
                                     igk_exx, evc_exx, transform_evc_to_exx
#if defined(__CUDA)
    USE device_memcpy_m,      ONLY : dev_memset
    USE device_fbuff_m,       ONLY : dev_buf
#endif
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ik, ibnd, i, j, k, ir, isym, ikq, ig, ierr
    INTEGER :: ibnd_loop_start
    INTEGER :: ipol, jpol
    REAL(DP), ALLOCATABLE :: occ(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: temppsic
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) temppsic
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    COMPLEX(DP),POINTER     :: psic_nc_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE)      :: psic_nc_d
#endif
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:), psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ibnd_start_new, ibnd_end_new, max_buff_bands_per_egrp
    INTEGER :: ibnd_exx, evc_offset
    !
    CALL start_clock ('exxinit2')
    !
    !$acc update device(evc)
    CALL transform_evc_to_exx( 2 )
    !
    ! Note that nxxs is not the same as nrxxs in parallel case
    nxxs = dfftt%nr1x * dfftt%nr2x * dfftt%nr3x
    nrxxs = dfftt%nnr
    !
#if defined(__MPI)
    IF (noncolin) THEN
       ALLOCATE( psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( temppsic(nrxxs) )
    ENDIF
    !
    CALL divide( inter_egrp_comm, x_nbnd_occ, ibnd_start, ibnd_end )
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    ! ... this will cause exxbuff to be calculated for every band
    ibnd_start_new = iexx_start
    ibnd_end_new = iexx_end
    !
    IF ( gamma_only ) THEN
        ibnd_buff_start = (ibnd_start_new+1)/2
        ibnd_buff_end   = (ibnd_end_new+1)/2
        max_buff_bands_per_egrp = MAXVAL((all_end(:)+1)/2-(all_start(:)+1)/2)+1
    ELSE
        ibnd_buff_start = ibnd_start_new
        ibnd_buff_end   = ibnd_end_new
        max_buff_bands_per_egrp = MAXVAL(all_end(:)-all_start(:))+1
    ENDIF
    !
    IF (.NOT. ALLOCATED(exxbuff)) THEN
       IF (gamma_only) THEN
          ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                        max_buff_bands_per_egrp-1,nkqs) ) ! THIS WORKS as for k
       ELSE
          ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                        max_buff_bands_per_egrp-1,nkqs) )
       ENDIF
    ENDIF
    !
    IF (.not. allocated(exxbuff_d) .and. use_gpu) THEN
       IF (gamma_only) THEN
          ALLOCATE( exxbuff_d(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nks))
       ELSE
          ALLOCATE( exxbuff_d(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nkqs))
       END IF
    ENDIF
    !
    IF (use_gpu) THEN
#if defined (__CUDA)
       ! NB: the array bounds are not passed to the subroutine.
       !
       ! See https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/269311
       !
       ! NB: TO BE CORRECTED WITH THE NEW DeviceXlib LIBRARY that dues internal slicing right!
       CALL dev_memset(exxbuff_d, (0.0_DP,0.0_DP), &
                                 (/ 1,nrxxs*npol/), 1, &
                                 (/ ibnd_buff_start, ibnd_buff_end /), ibnd_buff_start, &
                                 (/ 1,SIZE(exxbuff_d,3)/), 1)
#endif
    ELSE
       !$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs, &
       !$omp                ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
       DO ikq = 1, SIZE(exxbuff,3) 
          DO ibnd = ibnd_buff_start, ibnd_buff_end
             DO ir = 1, nrxxs*npol
                exxbuff(ir,ibnd,ikq) = (0.0_DP,0.0_DP)
             ENDDO
          ENDDO
       ENDDO
       ! the above loops will replaced with the following line soon
       !CALL threaded_memset(exxbuff, 0.0_DP, nrxxs*npol*SIZE(exxbuff,2)*nkqs*2)
    ENDIF
    !
    ! ... This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks
       !
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          IF (MOD(iexx_start,2) == 0) THEN
             ibnd_loop_start = iexx_start-1
          ELSE
             ibnd_loop_start = iexx_start
          ENDIF
          !
          evc_offset = 0
          DO ibnd = ibnd_loop_start, iexx_end, 2
             !
             psic(:) = ( 0._DP, 0._DP )
             !
             IF ( ibnd < iexx_end ) THEN
                IF ( ibnd == ibnd_loop_start .AND. MOD(iexx_start,2) == 0 ) THEN
                   DO ig = 1, npwt
                      psic(dfftt%nl(ig))  = ( 0._DP, 1._DP )*evc_exx(ig,1)
                      psic(dfftt%nlm(ig)) = ( 0._DP, 1._DP )*CONJG(evc_exx(ig,1))
                   ENDDO
                   evc_offset = -1
                ELSE
                   DO ig = 1, npwt
                      psic(dfftt%nl(ig))  = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) &
                           + ( 0._DP, 1._DP ) * evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2)
                      psic(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) ) &
                           + ( 0._DP, 1._DP ) * CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2) )
                   ENDDO
                ENDIF
             ELSE
                DO ig=1,npwt
                   psic(dfftt%nl (ig)) = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)
                   psic(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) )
                ENDDO
             ENDIF
             !
             CALL invfft( 'Wave', psic, dfftt )
             !
             exxbuff(1:nrxxs,(ibnd+1)/2,current_ik)=psic(1:nrxxs) 
             !
          ENDDO
          !
       ELSE IF_GAMMA_ONLY
          !
          npw = ngk (ik)
          IBND_LOOP_K : &
          DO ibnd = iexx_start, iexx_end
             !
             ibnd_exx = ibnd
             IF (noncolin) THEN
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir = 1, nrxxs
                   temppsic_nc(ir,1) = ( 0._DP, 0._DP )
                   temppsic_nc(ir,2) = ( 0._DP, 0._DP )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig = 1, npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),1) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic_nc(:,1), dfftt )
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx,npwx)
                DO ig = 1, npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),2) = evc_exx(ig+npwx,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic_nc(:,2), dfftt )
             ELSE
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir = 1, nrxxs
                   temppsic(ir) = ( 0._DP, 0._DP )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig = 1, npw
                   temppsic(dfftt%nl(igk_exx(ig,ik))) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic, dfftt )
             ENDIF
             !
             DO ikq = 1, nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = ABS(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#if defined(__MPI)
                   DO ipol = 1, npol
                      CALL gather_grid( dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol) )
                   ENDDO
                   !
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do collapse(2)
                      DO ipol = 1, npol
                         DO ir = 1, nxxs
                            psic_all_nc(ir,ipol) = (0.0_DP, 0.0_DP)
                            DO jpol = 1, npol
                               psic_all_nc(ir,ipol) = psic_all_nc(ir,ipol) + &
                                             CONJG(d_spin(jpol,ipol,isym)) * &
                                             temppsic_all_nc(rir(ir,isym),jpol)
                            ENDDO
                         ENDDO
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   DO ipol = 1, npol
                      CALL scatter_grid( dfftt, psic_all_nc(:,ipol), psic_nc(:,ipol) )
                   ENDDO
#else
!$omp parallel do collapse(2)
                   DO ipol = 1, npol
                      DO ir = 1, nxxs
                         psic_nc(ir,ipol) = (0._DP,0._DP)
                         DO jpol = 1, npol
                            psic_nc(ir,ipol) = psic_nc(ir,ipol) + CONJG(d_spin(jpol,ipol,isym))* &
                                               temppsic_nc(rir(ir,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDDO
!$omp end parallel do
#endif
                   !
#if defined (__CUDA)
                   IF (use_gpu) CALL dev_buf%lock_buffer(psic_nc_d, (/nrxxs, npol/), ierr)
                   IF (use_gpu) psic_nc_d = psic_nc
#endif
                   !
                   IF (index_sym(ikq) > 0 ) THEN
                      IF (use_gpu) THEN
                         associate(exxbuff=>exxbuff_d, psic_nc=>psic_nc_d)
                         ! sym. op. without time reversal: normal case
                         !$cuf kernel do 
                         DO ir=1,nrxxs
                            exxbuff(ir,ibnd,ikq)=psic_nc(ir,1)
                            exxbuff(ir+nrxxs,ibnd,ikq)=psic_nc(ir,2)
                         ENDDO
                         end associate
                      ELSE
                         ! sym. op. without time reversal: normal case
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir = 1, nrxxs
                         exxbuff(ir,ibnd,ikq) = psic_nc(ir,1)
                         exxbuff(ir+nrxxs,ibnd,ikq) = psic_nc(ir,2)
                      ENDDO
!$omp end parallel do
                      END IF
                   ELSE
                      ! sym. op. with time reversal: spin 1->2*, 2->-1*
                      IF (use_gpu) THEN
                         associate(exxbuff=>exxbuff_d, psic_nc=>psic_nc_d)
                         ! sym. op. with time reversal: spin 1->2*, 2->-1*
                         !$cuf kernel do 
                         DO ir=1,nrxxs
                            exxbuff(ir,ibnd,ikq)=CONJG(psic_nc(ir,2))
                            exxbuff(ir+nrxxs,ibnd,ikq)=-CONJG(psic_nc(ir,1))
                         ENDDO
                         end associate
                      ELSE
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir = 1, nrxxs
                         exxbuff(ir,ibnd,ikq) = CONJG(psic_nc(ir,2))
                         exxbuff(ir+nrxxs,ibnd,ikq) = -CONJG(psic_nc(ir,1))
                      ENDDO
!$omp end parallel do
                      ENDIF
                   ENDIF
#if defined(__CUDA)
                IF (use_gpu) CALL dev_buf%release_buffer(psic_nc_d, ierr)
                IF (use_gpu) exxbuff = exxbuff_d
#endif
                ELSE ! noncolinear
#if defined(__MPI)
                   CALL gather_grid( dfftt, temppsic, temppsic_all )
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                      DO ir = 1, nxxs
                         psic_all(ir) = temppsic_all(rir(ir,isym))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   CALL scatter_grid( dfftt, psic_all, psic )
#else
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                   DO ir = 1, nrxxs
                      psic(ir) = temppsic(rir(ir,isym))
                   ENDDO
!$omp end parallel do
#endif
!$omp parallel do default(shared) private(ir) firstprivate(isym,ibnd,ikq)
                   DO ir = 1, nrxxs
                      IF (index_sym(ikq) < 0 ) THEN
                         psic(ir) = CONJG(psic(ir))
                      ENDIF
                      exxbuff(ir,ibnd,ikq) = psic(ir)
                   ENDDO
!$omp end parallel do
                   !
                ENDIF ! noncolinear
                !
             ENDDO
             !
          ENDDO&
          IBND_LOOP_K
          !
       ENDIF&
       IF_GAMMA_ONLY
    ENDDO&
    KPOINTS_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, psic_nc )
#if defined(__MPI)
       DEALLOCATE( temppsic_all_nc, psic_all_nc )
#endif
    ELSE IF ( .NOT. gamma_only ) THEN
       DEALLOCATE( temppsic )
#if defined(__MPI)
       DEALLOCATE( temppsic_all, psic_all )
#endif
    ENDIF
    !
    ! Each wavefunction in exxbuff is computed by a single pool, collect among 
    ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
    ! See also the initialization of working_pool in exx_mp_init
    ! Note that in Gamma-only LSDA can be parallelized over two pools, and there
    ! is no need to communicate anything: each pools deals with its own spin
    !
    IF ( .NOT. gamma_only ) THEN
       DO ikq = 1, nkqs
         CALL mp_bcast( exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm ) 
       ENDDO
    ENDIF
    !
    ! For US/PAW only: compute <beta_I|psi_j,k+q> for the entire 
    ! de-symmetrized k+q grid by rotating the ones from the irreducible wedge
    !
    IF (okvan) CALL rotate_becxx( nkqs, index_xk, index_sym, xkq_collect )
    !
    ! Initialize 4-wavefunctions one-center Fock integrals
    !    \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
    !
    IF (okpaw) CALL PAW_init_fock_kernel()
    !
    CALL change_data_structure( .FALSE. )
    !
    CALL stop_clock( 'exxinit2' )
    !
  END SUBROUTINE exxinit2
  !
END MODULE exx2
