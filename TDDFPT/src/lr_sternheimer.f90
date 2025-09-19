!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE lr_sternheimer
  !
  !    This routine generalizes to finite complex frequencies and
  !    finite q vectors the routine solve_e of the Quantum ESPRESSO
  !    distribution.
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field
  !    of finite wavevector q and complex frequency omega.
  !    It performs the following tasks:
  !     a) computes the bare potential term  e^{iqr} | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system at zero
  !        frequency or ccg_many_vectors
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !
CONTAINS

SUBROUTINE one_sternheimer_step(iu, flag)
    !
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : e2, fpi, rytoev
    USE ions_base,              ONLY : nat
    USE io_global,              ONLY : stdout, ionode
    USE io_files,               ONLY : diropn, nwordwfc
    USE cell_base,              ONLY : tpiba2
    USE fft_interfaces,         ONLY : fwfft
    USE fft_interfaces,         ONLY : fft_interpolate
    USE klist,                  ONLY : lgauss, xk, wk
    USE gvecs,                  ONLY : doublegrid
    USE fft_base,               ONLY : dfftp, dffts
    USE lsda_mod,               ONLY : lsda, nspin, current_spin, isk
    USE wvfct,                  ONLY : nbnd, npwx, et
    USE klist,                  ONLY : ngk, igk_k
    USE check_stop,             ONLY : check_stop_now
    USE buffers,                ONLY : get_buffer, save_buffer
    USE wavefunctions,          ONLY : evc
    USE uspp,                   ONLY : vkb
    USE uspp_param,             ONLY : nhm
    USE noncollin_module,       ONLY : noncolin, domag, npol, nspin_mag
    USE scf,                    ONLY : rho
    USE gvect,                  ONLY : gg
    USE paw_variables,          ONLY : okpaw
    USE paw_onecenter,          ONLY : paw_dpotential
    USE eqv,                    ONLY : dpsi, dvpsi, evq
    USE units_lr,               ONLY : lrwfc, iuwfc, iudwf, lrdwf
    USE control_lr,             ONLY : lgamma, alpha_pv, nbnd_occ, &
                                       ext_recover, rec_code, &
                                       lnoloc, convt, tr2_ph, &
                                       alpha_mix, lgamma_gamma, niter_ph, &
                                       flmixdpot, rec_code_read

    USE lrus,                   ONLY : int3_paw
    USE qpoint,                 ONLY : xq, nksq, ikks, ikqs
    USE linear_solvers,         ONLY : ccg_many_vectors
    USE dv_of_drho_lr,          ONLY : dv_of_drho
    USE mp_bands,               ONLY : intra_bgrp_comm
    USE mp_images,              ONLY : root_image, my_image_id
    USE mp,                     ONLY : mp_sum
    USE fft_helper_subroutines, ONLY : fftx_ntgrp
    USE lr_variables,           ONLY : fru, fiu, iundvpsi, &
                                       iudrho, n_ipol, lr_verbosity, &
                                       chirr, chirz, chizr, chizz, epsm1, &
                                       current_w, itermax
    USE paw_add_symmetry,       ONLY : paw_deqsymmetrize
    USE wavefunctions,          ONLY : psic
    USE apply_dpot_mod,         ONLY : apply_dpot_allocate, apply_dpot_deallocate, &
                                       apply_dpot_bands
    USE response_kernels,      ONLY : sternheimer_kernel, sternheimer_kernel_freq, &
                                      sternheimer_postprocess
    USE uspp_init,             ONLY : init_us_2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(IN) :: flag   ! if 1 compute the charge-charge and
                                  ! charge magnetization responses
                                  ! if 2 and lsda computes the magnetization
                                  ! magnetization response
    REAL(DP) ::  thresh, anorm, avg_iter, dr2
    ! thresh: convergence threshold
    ! anorm : the norm of the error
    ! avg_iter: average number of iterations
    ! dr2   : self-consistency error
    COMPLEX(DP), ALLOCATABLE :: h_diag (:,:)
    COMPLEX(DP), ALLOCATABLE :: h_diag1 (:,:)
    !
    COMPLEX(DP) , ALLOCATABLE ::   &
                   dpsi1(:,:),   &
                   drhoscfout (:,:), & ! change of the scf charge (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   aux2(:,:), dvpsi1(:,:)
    !
    LOGICAL :: conv_root
    ! conv_root: true if linear system is converged
    INTEGER :: kter, iter0, ibnd, iter, lter, ik, ikk, ikq, &
               is, nrec, ndim, npw, npwq
    ! counters
    INTEGER :: nsolv
    !! Number of linear systems to solve. (1 for zero frequency, 2 for finite frequency)
    INTEGER :: isolv
    !! Index of the linear system to solve.
    INTEGER :: tot_num_iter, tot_cg_calls
    REAL(DP) :: xqmod2, alpha_pv0
    !
    REAL(DP) :: tcpu, get_clock
    ! timing variables
    !
    COMPLEX(DP) :: w  !frequency
    LOGICAL :: ldpsi1
    !
    EXTERNAL ch_psi_all, cg_psi
    EXTERNAL ch_psi_all_complex, ccg_psi
    !
    ! Input variables in dfpt_kernels
    COMPLEX(DP), ALLOCATABLE :: drhos(:, :, :)
    !! change of the charge density (smooth part only, but allocated with dfftp)
    COMPLEX(DP), ALLOCATABLE :: drhop(:, :, :)
    !! change of the charge density (smooth and hard parts, dfftp)
    COMPLEX(DP), POINTER :: dvscfs(:, :, :)
    !! change of the scf potential (smooth part only, dffts)
    COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfp(:, :, :)
    !! change of the scf potential (smooth and hard parts, dfftp)
    !
    ! Local variables in dfpt_kernels
    INTEGER :: ndim_pot, ndim_paw
    COMPLEX(DP), ALLOCATABLE :: dvscftmp(:, :, :)
    !! change of scf potential (output before mixing)
    !
    ! Input variables in sternheimer_kernel
    LOGICAL :: time_reversed
    !! True if solving a time-reversed Sternheimer equation
    LOGICAL :: all_conv
    !! True if sternheimer_kernel is converged at all k points and perturbations
    LOGICAL :: first_iter
    !! true if the first iteration.
    INTEGER :: npert
    !! number of perturbations
    COMPLEX(DP), ALLOCATABLE :: drhoout(:, :, :)
    !! induced charge density
    !
    ! Local variables in sternheimer_kernel
    INTEGER :: ipert
    REAL(DP) :: rsign
   !! sign of the term in the magnetization
    !
    CALL start_clock ('stern_step')
    !
    ! NOTE: In EELS, n_ipol = 1. In this code, n_ipol = 1 is assumed in a few places.
    !
    npert = n_ipol
    !
    w=CMPLX(fru(iu),fiu(iu))
    ldpsi1=ABS(w)>1.D-7
    alpha_pv0=alpha_pv
    alpha_pv=alpha_pv0 + REAL(w)
    !
    ALLOCATE (drhoscfout(dfftp%nnr, nspin_mag))
    !
    ndim_pot = dfftp%nnr * nspin_mag * 1
    IF (okpaw) THEN
       ndim_paw = (nhm * (nhm+1) * nat * nspin_mag * 1) / 2
       ALLOCATE(mixin(ndim_pot + ndim_paw))
       ALLOCATE(mixout(ndim_pot + ndim_paw))
       mixin = (0.0_DP, 0.0_DP)
    ELSE
       ALLOCATE(mixin(1))
       ALLOCATE(mixout(1))
    ENDIF
    !
    ALLOCATE (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 1))
    IF (noncolin) ALLOCATE (dbecsum_nc (nhm, nhm, nat, nspin, 1))
    IF (ldpsi1) THEN
       ALLOCATE (dpsi1(npwx*npol,nbnd))
       ALLOCATE (dvpsi1(npwx*npol,nbnd))
       ALLOCATE (h_diag(npwx*npol, nbnd))
       ALLOCATE (h_diag1(npwx*npol, nbnd))
    ENDIF
    ALLOCATE (aux2(npwx*npol, nbnd))
    ALLOCATE(dvscftmp(dfftp%nnr, nspin_mag, npert))
    ALLOCATE(drhos(dffts%nnr, nspin_mag, npert))
    ALLOCATE(drhop(dfftp%nnr, nspin_mag, npert))
    ALLOCATE(dvscfp(dfftp%nnr, nspin_mag, npert))
    dvscfp(:,:,:)=(0.d0,0.d0)
    IF (doublegrid) THEN
       ALLOCATE(dvscfs(dffts%nnr, nspin_mag, npert))
    ELSE
       dvscfs => dvscfp
    ENDIF
    ALLOCATE(drhoout(dffts%nnr, nspin_mag, npert))
    CALL apply_dpot_allocate()
    !
    !$acc enter data create(aux2, dvscfs)
    dvpsi =(0.0d0, 0.0d0)

!    IF (rec_code_read == -20.AND.ext_recover) then
!       ! restarting in Electric field calculation
!       IF (okpaw) THEN
!          CALL read_rec(dr2, iter0, 1, dvscfp, dvscfs, drhop, dbecsum)
!          CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
!                      mixin, dvscfp, dbecsum, ndim, -1 )
!       ELSE
!          CALL read_rec(dr2, iter0, 1, dvscfp, dvscfs)
!       ENDIF
!    ELSEIF (rec_code_read > -20 .AND. rec_code_read <= -10) then
!       ! restarting in Raman: proceed
!       convt = .true.
!    ELSE
       convt = .false.
       iter0 = 0
!    ENDIF
    !
!    IF ( ionode .AND. fildrho /= ' ') THEN
!       INQUIRE (UNIT = iudrho, OPENED = exst)
!       IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
!       CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
!    ENDIF
    IF (rec_code_read > -20) convt=.TRUE.
    !
    IF (convt) go to 155
    !
    IF ((lgauss.and..not.ldpsi1)) &
            CALL errore ('solve_eq', 'insert a finite frequency', 1)
    !
    IF (lr_verbosity > 5) THEN
       WRITE(stdout,'("<lr_sternheimer_one_step>")')
    ENDIF
    !
    IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
    !
    IF (ldpsi1) THEN
       nsolv = 2  ! Finite frequency, solve original and time-reversed Sternheimer equation
    ELSE
       nsolv = 1  ! Zero frequency, solve only one Sternheimer equation
       ! For noncollinear magnetism, one needs to solve two Sternheimer equations even for
       ! zero frequency. This is not implemented here.
    ENDIF
    !
    ! Calculate bare perturbation multiplied to the wavefunctions, save on buffer iundvpsi
    !
    DO ik = 1, nksq
       !
       ikk  = ikks(ik)
       ikq  = ikqs(ik)
       npw  = ngk(ikk)
       npwq = ngk(ikq)
       IF (lsda) current_spin = isk (ikk)
       !
       ! Read unperturbed wavefuctions evc (wfct at k)
       ! and evq (wfct at k+q)
       !
       IF (nksq > 1) THEN
          CALL get_buffer(evc, lrwfc, iuwfc, ikk)
       ENDIF
       !
       ! Calculate beta-functions vkb at k+q (Kleinman-Bylander projectors)
       ! The vkb's are needed for the non-local potential in h_psi,
       ! and for the ultrasoft term.
       !
       CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
       !$acc update host(vkb)
       !
       ! do over polarization
       !
       DO ipert = 1, npert
          !
          nrec = (ipert - 1) * nksq + ik
          !  IF (isolv == 2) nrec = nrec + npert * nksq
          !
          CALL dveqpsi_us(ik)
          !
          !  with flag=2 the perturbation is a magnetic field along z
          !
          IF (lsda.AND.current_spin==2.AND.flag==2) dvpsi=-dvpsi
          !
          CALL save_buffer(dvpsi, nwordwfc, iundvpsi, nrec)
          !
       ENDDO ! ipert
    ENDDO ! ik
    !
    ! Loop over the iterations
    !
    DO kter = 1, itermax
       !
       FLUSH( stdout)
       iter = kter + iter0
       first_iter = (iter == 1)
       !
       tot_num_iter = 0
       tot_cg_calls = 0
       !
       drhos(:,:,:) = (0.0d0, 0.0d0)
       dbecsum(:,:,:,:) = (0.0d0, 0.0d0)
       IF (noncolin) dbecsum_nc = (0.0d0, 0.0d0)
       !
       ! Set threshold for iterative solution of the linear system (ccgsolve_all)
       !
       IF (first_iter) THEN
          thresh = 1.d-2
          IF (lnoloc) thresh = 1.d-5
       ELSE
          thresh = min(1.d-1 * sqrt(dr2), thresh)
       ENDIF
       !
       ! Compute drhos, the charge density response to the total potential
       !
       time_reversed = .FALSE.  ! magnetic case not implemented
       !
       IF (ldpsi1) THEN
          !
          ! Finite-frequency Sternheimer equation
          CALL sternheimer_kernel_freq(first_iter, time_reversed, npert, nwordwfc, iundvpsi, &
             thresh, w, dvscfs, all_conv, avg_iter, drhos, dbecsum, dbecsum_nc)
          !
       ELSE
          !
          ! Zero-frequency Sternheimer equation
          CALL sternheimer_kernel(first_iter, time_reversed, npert, nwordwfc, iundvpsi, &
             thresh, dvscfs, all_conv, avg_iter, drhos, dbecsum, dbecsum_nc)
          !
       ENDIF
       !
       ! Postprocess the results of the Sternheimer kernel
       !
       CALL sternheimer_postprocess(nsolv, npert, drhos, drhop, dbecsum, dbecsum_nc)
       !
       !   compute the corresponding change in scf potential : drhop -> dvscftmp
       !
       IF (lnoloc) THEN
         ! No local field effect: set dvscf to 0
         dvscftmp(:, :, ipert) = (0.d0, 0.d0)
      ELSE
         ! Compute the response HXC potential
         CALL zcopy(dfftp%nnr*nspin_mag, drhop(1,1,1), 1, dvscftmp(1,1,1), 1)
         CALL dv_of_drho(dvscftmp(1, 1, 1))
       ENDIF
       !
       !  mix with the old potential: dvscftmp -> dvscfp
       !
       IF (okpaw) THEN
          !
          !  In this case we mix also dbecsum
          !
          CALL setmixout(ndim_pot, ndim_paw, mixout, dvscftmp, dbecsum, ndim, -1)
          CALL mix_potential_eels(2*(ndim_pot + ndim), mixout, mixin, alpha_mix(kter), dr2, &
                                  tr2_ph / npol, iter, flmixdpot, convt)
          CALL setmixout(ndim_pot, ndim_paw, mixin, dvscfp, dbecsum, ndim, 1)
       ELSE
          ! nmix_ph ??
          CALL mix_potential_eels(2*ndim_pot, dvscftmp, dvscfp, alpha_mix(kter), dr2, &
                                  tr2_ph / npol, iter, flmixdpot, convt)
       ENDIF
       !
       IF (doublegrid) then
          DO is = 1, nspin_mag
             CALL fft_interpolate(dfftp, dvscfp(:, is, 1), dffts, dvscfs(:, is, 1))
          ENDDO
       ENDIF
       !
       IF (okpaw) THEN
          IF (noncolin.AND.domag) THEN
!             CALL PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
          ELSE
             !
             ! The presence of c.c. in the formula gives a factor 2.0
             !
             dbecsum=2.0_DP * dbecsum
             IF (.NOT. lgamma_gamma) CALL paw_deqsymmetrize(dbecsum)
             CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,1)
          ENDIF
       ENDIF
       !
       CALL newdq(dvscfp, 1)
       !
!       tcpu = get_clock ('PHONON')
       tcpu = get_clock ('ccgsolve')
       WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
            &      " secs   av.it.: ",f5.1)') iter, tcpu, avg_iter
       WRITE( stdout, "(5x,' thresh=',es10.3, ' alpha_mix = ',f6.3, &
            &      ' |ddv_scf|^2 = ',es10.3 )") thresh, alpha_mix (kter), dr2
       !
       FLUSH( stdout )
       !
       ! rec_code: state of the calculation
       ! rec_code=-20 Electric Field
       !
       rec_code=-20
!       IF (okpaw) THEN
!          CALL write_rec('solve_e...', 0, dr2, iter, convt, 1, dvscfp, &
!                                                      drhop, dbecsum)
!       ELSE
!          CALL write_rec('solve_e...', 0, dr2, iter, convt, 1, dvscfp)
!       ENDIF
       !
!       IF (check_stop_now()) CALL stop_smoothly_ph (.false.)
       !
       IF (convt) goto 155
       !
    ENDDO  ! do over iter_ph
    !
155 CONTINUE
    !
    drhoscfout(:,:) = drhop(:,:,1)
    !
    !  compute here the susceptibility and the inverse of the dielectric
    !  constant
    !
    !  CALL compute_susceptibility(drhoscfout)
    !
    DO is=1,nspin_mag
       CALL fwfft ('Rho', drhoscfout(:,is), dfftp)
    ENDDO
    !
    IF (flag==1) THEN
       chirr(iu)=(0.0_DP,0.0_DP)
       chizr(iu)=(0.0_DP,0.0_DP)
       epsm1(iu)=(0.0_DP,0.0_DP)
    ELSE
       chirz(iu)=(0.0_DP,0.0_DP)
       chizz(iu)=(0.0_DP,0.0_DP)
    ENDIF
    !
    xqmod2=(xq(1)**2+xq(2)**2+xq(3)**2)*tpiba2
    !
    IF (ABS(gg(1))<1.d-8) THEN
       IF (flag==1) THEN
          chirr(iu) = drhoscfout(dfftp%nl(1),1)
          IF (lsda) chirr(iu) = chirr(iu) + drhoscfout(dfftp%nl(1),2)
          epsm1(iu) = CMPLX(1.0_DP,0.0_DP)+ chirr(iu)*fpi*e2/xqmod2
          IF (lsda) chizr(iu) = drhoscfout(dfftp%nl(1),1) - &
                                drhoscfout(dfftp%nl(1),2)
       ELSEIF (lsda) THEN
          chizz(iu)=drhoscfout(dfftp%nl(1),1)-drhoscfout(dfftp%nl(1),2)
          chirz(iu)=drhoscfout(dfftp%nl(1),1)+drhoscfout(dfftp%nl(1),2)
       ENDIF
    ENDIF
    !
    IF (flag==1) THEN
       CALL mp_sum(epsm1(iu),intra_bgrp_comm)
       CALL mp_sum(chirr(iu),intra_bgrp_comm)
       CALL mp_sum(chizr(iu),intra_bgrp_comm)
    ELSE
       CALL mp_sum(chizz(iu),intra_bgrp_comm)
       CALL mp_sum(chirz(iu),intra_bgrp_comm)
    ENDIF
    !
    IF (flag==1) THEN
       WRITE(stdout, '(/,6x,"Inverse dielectric constant at &
                          &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
       WRITE(stdout, '(46x,f9.4," +",f9.4," i eV")') current_w * rytoev
       WRITE(stdout,'(/,6x,"epsilon^-1(q,w) =",2f15.6)') epsm1(iu)
       !
       WRITE( stdout, '(/,6x,"Charge-charge susceptibility:")')
       !
       WRITE(stdout,'(/,6x,"chirr(q,w) =",2f15.6)') chirr(iu)
       IF (lsda) THEN
          WRITE(stdout,'(/,6x,"m_z-charge susceptibility:")')
          WRITE(stdout,'(/,6x,"chizr(q,w) =",2f15.6)') chizr(iu)
       ENDIF
       !
    ELSEIF (lsda) THEN
       WRITE( stdout, '(/,6x,"m_z - m_z susceptibility at &
                       &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
       WRITE( stdout, '(43x,f9.4," +",f9.4," i eV")') current_w * rytoev
       WRITE(stdout,'(/,6x,"chizz(q,w) =",2f15.6)') chizz(iu)
       WRITE(stdout,'(/,6x,"chirz(q,w) =",2f15.6)') chirz(iu)
    ENDIF
    !
    CALL apply_dpot_deallocate()
    IF (ldpsi1) THEN
       deallocate (dpsi1)
       deallocate (dvpsi1)
       deallocate (h_diag)
       deallocate (h_diag1)
    ENDIF
    deallocate (dbecsum)
    IF (okpaw) THEN
       DEALLOCATE(mixin)
       DEALLOCATE(mixout)
    ENDIF
    deallocate (drhoscfout)
    !$acc exit data delete(aux2, dvscfs)
    if (noncolin) deallocate(dbecsum_nc)
    deallocate(aux2)
    IF (doublegrid) DEALLOCATE (dvscfs)
    DEALLOCATE(drhos)
    DEALLOCATE(drhop)
    DEALLOCATE(dvscfp)
    DEALLOCATE(drhoout)
    !
    alpha_pv=alpha_pv0
    !
    CALL stop_clock ('stern_step')
    !
    RETURN
    !
END SUBROUTINE one_sternheimer_step

END MODULE lr_sternheimer
