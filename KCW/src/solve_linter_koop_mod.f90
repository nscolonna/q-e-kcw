! Copyright (C) 2003-2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE solve_linter_koop_mod

CONTAINS
!
!!JA ASSUME all arrays are on GPU excpet
!!JA drhog_scf and drhor_scf - are crateted and return ONLY on GPU
!!JA delta_vg is only used for OLD  - sternheimer_kernel_old - not USED
!-----------------------------------------------------------------------
subroutine solve_linter_koop ( spin_ref, i_ref, delta_vr, drhog_scf, delta_vg, drhor_scf)
  !-----------------------------------------------------------------------
  !
  !!    Driver routine for the solution of the linear system which defines the 
  !!    change of the wavefunction due to an EXTERNAL perturbation. It work 
  !!    exactly as solve_linter in PH but for a genenral perturbation dv delta_vr.
  !!    In general, it performs the following tasks:
  !!     a) computes the bare potential term Delta V | psi >, add the scf contribution
  !!     b) applies P_c^+ (orthogonalization to valence states)
  !!     c) calls cgstabsolve_all to solve the linear system
  !!     d) computes Delta rho
  !
  ! ## FIXME delta_vg passed only for debug reason (to test symmetries) Needs to be removed once the
  !          problem will be solved
  !
  USE control_kcw
  USE dv_of_drho_lr
  USE kinds,                 ONLY : DP
  USE constants,             ONLY : tpi
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout
  USE wavefunctions,         ONLY : evc, psic
  USE klist,                 ONLY : lgauss, ngk
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE fft_base,              ONLY : dffts, dfftp
  USE fft_interfaces,        ONLY : fwfft, invfft, fft_interpolate

  USE gvect,                 ONLY : gstart
  USE gvecs,                 ONLY : doublegrid, ngms
  USE becmod,                ONLY : calbec
  USE wvfct,                 ONLY : npw,npwx
  USE uspp_param,            ONLY : nhm
  USE control_lr,            ONLY : lgamma
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE buffers,               ONLY : save_buffer, get_buffer
  USE eqv,                   ONLY : dvpsi
  USE qpoint,                ONLY : nksq, ikks
  USE uspp,                  ONLY : okvan  
  USE uspp_init,             ONLY : init_us_2
  USE mp,                    ONLY : mp_sum
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_mag, npol
  USE mp_pools,              ONLY : inter_pool_comm
  USE apply_DPot_mod,        ONLY : apply_DPot_bands
  USE response_kernels,      ONLY : sternheimer_kernel
  USE scf,                   ONLY : vrs
  USE qpoint_aux,            ONLY : ikmks
  USE lr_symm_base,          ONLY : lr_npert, upert, upert_mq, minus_q, nsymq
  USE qpoint,                ONLY : xq
  USE symm_base,             ONLY : ft
  USE cell_base,             ONLY : at
  ! For LDOS at q=0 - apparently unused
  USE dfpt_type,             ONLY : dfpt_ldos_type, allocate_dfpt_ldos, &
                                    deallocate_dfpt_ldos
  USE localdos_mod,          ONLY : localdos_new
  USE ener,                  ONLY : ef
  ! end LDOS at q=0 
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout) ::delta_vr (dffts%nnr, nspin_mag)
  COMPLEX(DP), INTENT(out) ::  drhog_scf (ngms, nspin_mag)
  COMPLEX(DP), INTENT(out), optional :: drhor_scf(dffts%nnr,nspin_mag)
  !
  ! For LDOS at q=0 - apparently unused
  TYPE(dfpt_ldos_type) :: ldos_data
  ! end LDOS at q=0 
  REAL(DP) :: averlt, &
              thresh, & ! convergence threshold
              dr2       ! self-consistency error
  !
  COMPLEX(DP), ALLOCATABLE :: & 
              drhoscf (:,:),    & !
              drhoscfh (:,:),   & !
              dvscfout (:,:),   & !
              dbecsum(:,:,:,:), & !
              dbecsum_nc (:,:,:,:,:,:), & !
              aux(:)           !
  INTEGER :: nrec
  !
  COMPLEX(DP), POINTER ::  dvscfin(:,:,:), dvscfins (:,:,:)
  ! change of the scf potential,  change of the scf potential (smooth part only)
  REAL(DP), ALLOCATABLE :: becsum1(:,:,:)
  !
  LOGICAL :: lmetq0,     & ! true if xq=(0,0,0) in a metal
             all_conv,   & ! true if the Linear system converged for all the bands and k-points
             convt         ! true if the SCF problem converged

  INTEGER :: iter,          & ! counter on iterations
             ik, ikk, ikmk, & ! counter on k points
             is,            & ! counter on spin polarizations
             spin_ref,      & ! the spin of the reference orbital
             i_ref,         & ! the orbital we want to keep fix
             isolv,         & ! counter on linear system
             nsolv           ! number of linear systems
  !
  REAL(DP) :: tcpu, get_clock ! timing variables
  EXTERNAL ch_psi_all,  cg_psi
  CHARACTER(LEN=256) :: flmixDPot = 'mixd'
  !
  COMPLEX(DP) :: delta_vg(ngms,nspin_mag)
  INTEGER:: norb
  !
  LOGICAL :: new
  COMPLEX(DP)    ::   imag
  REAL(DP)       ::   xq_cryst(3)
  INTEGER        ::   isym,nnr

  ! Set to false to revert to the previous implementation of the Linear solver (consistent with QE7.1)
  new = .true.
  !
  CALL start_clock ('solve_linter')
  !
  nsolv=1
  IF (noncolin.AND.domag) nsolv=2
  !
  !
  nnr = dfftp%nnr
  IF (doublegrid) THEN
     ALLOCATE (dvscfins (dffts%nnr , nspin_mag , 1))
    !$acc enter data create(dvscfins)
    !$acc kernels present(dvscfins)
    dvscfins(:,:,:) = (0.D0, 0.D0)
    !$acc end kernels
    !$acc update host(dvscfins)
  ELSE
    ALLOCATE (dvscfin (dfftp%nnr, nspin_mag, 1)) 
    !$acc enter data create(dvscfin(1:nnr, 1:nspin_mag, 1))
    !$acc kernels present(dvscfin)
    dvscfin(:,:,:) = (0.D0, 0.D0)
    !$acc end kernels
    !$acc update host(dvscfin)
  ENDIF
  !
  ALLOCATE (drhoscf  (dffts%nnr, nspin_mag)) 
  ALLOCATE (drhoscfh (dfftp%nnr, nspin_mag))
  ALLOCATE (dvscfout (dfftp%nnr, nspin_mag))    
  ALLOCATE (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag , 1)) 
  !
  ! Set symmetry representation in lr_symm_base
  !
  lr_npert = 1
  ALLOCATE(upert(lr_npert, lr_npert, nsymq))
  xq_cryst = xq
  imag = (0.D0, 1.D0)
  CALL cryst_to_cart(1,xq_cryst,at,-1)
  DO isym = 1, nsymq
     upert(1, 1, isym) = EXP(- imag * tpi * dot_product( xq_cryst(:), ft(:, isym) ) ) 
  ENDDO
  IF (minus_q) THEN
     ALLOCATE(upert_mq(lr_npert, lr_npert))
     upert_mq(1, 1) = (1.d0, 0.d0)
  ENDIF ! minus_q
  ! 
  IF (noncolin) allocate (dbecsum_nc (nhm,nhm, nat , nspin_mag , 1, nsolv))
  ALLOCATE (aux ( dffts%nnr ))    
  !
  drhoscf  = CMPLX(0.D0, 0.D0, kind =DP)
  drhoscfh = CMPLX(0.D0, 0.D0, kind =DP)
  !
  IF (new .AND. fix_orb) THEN 
     CALL infomsg('kcw_solve_linter','WARNING: fix_orb not supported anymore; set it to FALSE')
     fix_orb= .FALSE.
  ENDIF 
  IF (kcw_at_ks .AND. fix_orb) WRITE (stdout, '("FREEZING ORBITAL #", i4, 3x , "spin", i4)') i_ref, spin_ref
  !
  ! if q=0 for a metal: ALLOCATE and compute local DOS at Ef (unused?)
  ! 
  lmetq0 = lgauss .AND. lgamma
  IF (lmetq0) THEN
     CALL allocate_dfpt_ldos(ldos_data)
     CALL localdos_new(ldos_data, ef)
  ENDIF
  !
  DO ik = 1, nksq 
     !
     ikk = ikks(ik)
     npw = ngk(ikk)
     !
     IF (lsda) current_spin = isk (ikk)
     !
     !
     DO isolv = 1, nsolv
       IF (isolv == 1) THEN
          ikmk = ikks(ik)
       ELSE
           ikmk = ikmks(ik)
       ENDIF
       !
       ! read unperturbed wavefunctions psi(k) and psi(k+q)
       !
       IF (nksq .GT. 1 .OR. nsolv == 2) THEN 
           CALL get_buffer (evc, lrwfc, iuwfc, ikmk)
           !$acc update device(evc)
       END IF    
       nrec = (isolv-1) * nksq + ik
       !
       ! ... Compute dv_bare*psi (dvspi) and store it 
       !     This depends on isolv [if isolv=2 delta_vr(2:4) -> - delta_vr(2:4)]
       IF (isolv==2) THEN 
          !$acc kernels present(delta_vr)
          delta_vr(:,2:4) = - delta_vr(:,2:4)
          !$acc end kernels
       END IF   
       CALL kcw_dvqpsi (ik, delta_vr, isolv)
       CALL save_buffer (dvpsi, lrdvwfc, iudvwfc, nrec)
       !WRITE(*,'("NICOLA dvpsi", I5, 3X, 3(F20.18,2x))') nrec, REAL(CONJG(dvpsi(npwx+1:npwx+3,1))*dvpsi(npwx+1:npwx+3,1))
       ! Restore the sign of delta_vr
       IF (isolv==2) THEN 
          !$acc kernels present(delta_vr)
          delta_vr(:,2:4) = - delta_vr(:,2:4)
          !$acc end kernels
       END IF   
       !
       IF (okvan) THEN
          CALL errore('solve_linter_koop', 'USPP not implemented yet', 1)
       ENDIF
       !
     ENDDO ! isolve 
     !
  ENDDO ! ik
  !
  !   Loop is over the iterations
  !
  dr2=0.d0
  !$acc enter data create(drhoscf)
  DO iter = 1, niter
     !
     !$acc kernels present(drhoscf)
     drhoscf = (0.d0, 0.d0)
     !$acc end kernels
     dbecsum = (0.d0, 0.d0)
     !
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     !
     IF (doublegrid) THEN
        DO isolv = 1, nsolv
          !
          IF (iter == 1 ) THEN
             thresh = 1.d-6
          ELSE
            thresh = min (1.d-2 * sqrt (dr2), 1.d-6)
          ENDIF
          ! 
          !
          CALL sternheimer_kernel(iter==1, isolv==2, 1, lrdvwfc, iudvwfc, &
            thresh, dvscfins, all_conv, averlt, drhoscf, dbecsum, &
            dbecsum_nc(:,:,:,:,:,isolv))
          !
        ENDDO
     ELSE
        DO isolv = 1, nsolv
          !
          IF (iter == 1 ) THEN
            thresh = 1.d-6
          ELSE
            thresh = min (1.d-2 * sqrt (dr2), 1.d-6)
          ENDIF
          ! 
          !
          CALL sternheimer_kernel(iter==1, isolv==2, 1, lrdvwfc, iudvwfc, &
            thresh, dvscfin, all_conv, averlt, drhoscf, dbecsum, &
            dbecsum_nc(:,:,:,:,:,isolv))
          !
          !
        ENDDO
     END IF
     !
     IF (nsolv==2) THEN
        !$acc kernels present(drhoscf)
        drhoscf = drhoscf / 2.0_DP
        !$acc end kernels
        dbecsum = dbecsum / 2.0_DP
        IF (noncolin) dbecsum_nc = dbecsum_nc / 2.0_DP
     ENDIF
     !
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     IF (noncolin) THEN
        CALL mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        CALL mp_sum ( dbecsum, intra_bgrp_comm )
     ENDIF
     !
     !
     !$acc update host(drhoscf)
     IF (doublegrid) THEN
        DO is = 1, nspin_mag
           CALL fft_interpolate (dffts, drhoscf(:,is), dfftp, drhoscfh(:,is))
        ENDDO
     ELSE
        CALL zcopy (nspin_mag*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     ENDIF
     !
     !
     ! Symmetrization of the response charge density.
     !
     IF (irr_bz) CALL psymdvscf (drhoscfh, dfftp)
     !
     !
     !    Now we compute for all perturbations the total charge and potential
     !
     !CALL addusddens (drhoscfh, dbecsum, irr, imode0, 1, 0)
     !
     !   Reduce the delta rho across pools
     !
     CALL mp_sum (drhoscf, inter_pool_comm)
     CALL mp_sum (drhoscfh, inter_pool_comm)
     !$acc update device(drhoscf)
     !
     !   ... save them on disk and
     !   compute the corresponding change in scf potential
     !
     !! CALL zcopy (dfftp%nnr*nspin_mag,drhoscfh(1,1),1,dvscfout(1,1),1) !<--- WHY (1,1)???
     dvscfout(:,:) = drhoscfh(:,:)

     ! NB: always CALL with imode=0 to avoid CALL to addcore in dv_of_drho for 
     !     nlcc pseudo. The CALL is not needed since we are not moving atoms!!
     !
     CALL dv_of_drho (dvscfout)  
     !
     ! ... On output in dvscfin we have the mixed potential
     !
     !! TMP REMOVED $acc update host(dvscfin)
     CALL mix_potential (2*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &   
                         alpha_mix(iter), dr2, tr2/npol, iter, &
                         nmix, flmixDPot, convt)
     !$acc update device(dvscfin) 
     !WRITE(mpime+1000, '(1i5,es10.3,1l1,1i5)') my_pool_id, dr2, convt, iter
     !
     ! check that convergent have been reached on ALL processors in this image
     CALL check_all_convt(convt)
     !
     IF (doublegrid) THEN
        DO is = 1, nspin_mag
           CALL fft_interpolate (dffts, drhoscf(:,is), dfftp, drhoscfh(:,is))
        ENDDO
     ENDIF
     !
     tcpu = get_clock ('KCW')
     !
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f9.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (iter) , dr2
     !
     !    Here we save the information for recovering the run from this poin
     ! 
     FLUSH( stdout )
     !
     IF (convt) EXIT
     !
  ENDDO  ! loop over iteration

  IF (present(drhor_scf)) THEN 
      !$acc enter data create(drhor_scf)
      !$acc kernels present(drhor_scf, drhoscf)
      drhor_scf = drhoscf(:,:)
      !$acc end kernels
  END IF    
  !$acc exit data delete(drhoscf)
  !
  ! The density variation in G-space
  !
  !$acc enter data create(drhog_scf, aux)
  DO is = 1, nspin_mag
     aux(:) = drhoscfh(:,is)
     !$acc update device(aux)
     !$acc host_data use_device(aux)
     CALL fwfft ('Rho', aux, dffts)
     !$acc end host_data
     !$acc kernels present(drhog_scf, aux, dffts, dffts%nl) 
     drhog_scf(:,is) = aux(dffts%nl(:))
     !$acc end kernels
  ENDDO
  !$acc exit data delete(aux)
  !
  !
  ! The induced density in G space
  !
  if (lmetq0) CALL deallocate_dfpt_ldos(ldos_data)
  DEALLOCATE (aux)
  DEALLOCATE (dbecsum)
  IF (noncolin) DEALLOCATE (dbecsum_nc)
  DEALLOCATE (drhoscf )
  DEALLOCATE (dvscfout)
  DEALLOCATE (drhoscfh)
  IF (doublegrid) THEN 
     !$acc exit data delete(dvscfins)
     DEALLOCATE (dvscfins)
  ELSE   
     !$acc exit data delete(dvscfin)
     DEALLOCATE (dvscfin)
  END IF
  DEALLOCATE (upert)
  IF (minus_q) DEALLOCATE(upert_mq) 
  ! 
  !
  CALL stop_clock ('solve_linter')
  !
  RETURN
  !
END SUBROUTINE solve_linter_koop


!------------------------------------------------------------------
SUBROUTINE check_all_convt( convt )
  !---------------------------------------------------------------
  !! Work out how many processes have converged.
  !
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : nproc_image, intra_image_comm
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(in) :: convt
  INTEGER            :: tot_conv
  !
  IF(nproc_image==1) RETURN
  !
  tot_conv = 0
  IF(convt) tot_conv = 1
  CALL mp_sum(tot_conv, intra_image_comm)
  !
  IF ((tot_conv > 0) .and. (tot_conv < nproc_image)) THEN
    CALL errore('check_all_convt', 'Only some processors converged: '&
               &' either something is wrong with solve_linter, or a different'&
               &' parallelism scheme should be used.', 1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE


!!JA asume delta_vr are on GPU
!------------------------------------------------------------------
SUBROUTINE kcw_dvqpsi (ik, delta_vr, isolv)
  !----------------------------------------------------------------
  ! 
  ! Apply the bare perturbing potential to the unperturbed KS wfc. 
  !
  USE kinds,                 ONLY : DP
  USE eqv,                   ONLY : dvpsi
  USE control_lr,            ONLY : nbnd_occ
  USE wvfct,                 ONLY : npw, npwx
  USE fft_base,              ONLY : dffts
  USE klist,                 ONLY : igk_k, ngk
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE fft_wave,              ONLY : invfft_wave, fwfft_wave
  USE qpoint,                ONLY : npwq, ikks, ikqs
  USE lsda_mod,              ONLY : nspin, current_spin
  USE wavefunctions,         ONLY : evc
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_mag, npol
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik, isolv
  INTEGER :: ibnd, ig, ir, ikk, ikq, ip 
  COMPLEX(DP), INTENT(in) ::delta_vr (dffts%nnr, nspin_mag)
  !!COMPLEX(DP) :: aux( dffts%nnr,npol)!, aux_g(npwx*npol)
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:)  :: aux  !, aux_g(npwx*npol)
  COMPLEX(DP) :: sup, sdwn    
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  ALLOCATE (aux(dffts%nnr, npol))

  !$acc enter data create(dvpsi, aux) 
  !$acc kernels present(dvpsi)
    dvpsi(:,:) = (0.0_dp, 0.0_dp)
  !$acc end kernels
 !! aux(:,:) = (0.D0, 0.D0)
  !
  DO ibnd = 1, nbnd_occ (ik)
      !$acc kernels present(aux)
      aux(:,:) = (0.0_dp, 0.0_dp)
      !$acc end kernels

      IF (nspin==2 .OR. nspin==1) THEN
  !OLD FFT
  !       DO ig = 1, npw
  !          aux (dffts%nl(igk_k(ig,ikk)),1)=evc(ig,ibnd)
  !       ENDDO
  !       CALL invfft ('Wave', aux, dffts)
      !!   !$acc enter data copyin(aux)
         CALL invfft_wave (npwx, npw, igk_k (1,ikk), evc(:,ibnd), aux )
      !!   !$acc exit data copyout(aux)

         !$acc parallel loop present(aux, delta_vr)
         DO ir = 1, dffts%nnr
             aux(ir,1)=aux(ir,1)*delta_vr(ir,current_spin) 
         ENDDO
         !
 !        CALL fwfft ('Wave', aux, dffts)
 !        DO ig = 1, npwq
 !           dvpsi(ig,ibnd)=aux(dffts%nl(igk_k(ig,ikq)),1)
 !        ENDDO
    !!     !$acc enter data copyin(dvpsi, aux) 
         CALL fwfft_wave (npwx, npwq, igk_k (1,ikq), dvpsi(:,ibnd) , aux)
    !!     !$acc exit data copyout(dvpsi, aux)

      ELSEIF (nspin==4) THEN
       !!  !$acc enter data copyin(aux)
         CALL invfft_wave (npwx, npw, igk_k (1,ikk), evc(:,ibnd), aux )
       !!  !$acc exit data copyout(aux)

         IF (domag) then
            !$acc parallel loop present(aux, delta_vr) private(sup, sdwn)
            DO ir = 1, dffts%nnr
               sup=aux(ir,1)*(delta_vr(ir,1)+delta_vr(ir,4))+ &
                   aux(ir,2)*(delta_vr(ir,2)-(0.d0,1.d0)*delta_vr(ir,3))
               sdwn=aux(ir,2)*(delta_vr(ir,1)-delta_vr(ir,4)) + &
                    aux(ir,1)*(delta_vr(ir,2)+(0.d0,1.d0)*delta_vr(ir,3))
               aux(ir,1)=sup
               aux(ir,2)=sdwn
            ENDDO
         ELSE
            !$acc parallel loop present(aux, delta_vr) 
            DO ir = 1, dffts%nnr
               aux(ir,1)=aux(ir,1)*delta_vr(ir,1)
               aux(ir,2)=aux(ir,2)*delta_vr(ir,1)
            ENDDO
         END IF
      !!   !$acc enter data copyin(dvpsi, aux) 
         CALL fwfft_wave (npwx, npwq, igk_k (1,ikq),dvpsi(:,ibnd) , aux)
      !!   !$acc exit data copyout(dvpsi, aux)
      ENDIF
     !
  ENDDO
  !$acc exit data copyout(dvpsi) delete(aux)
  DEALLOCATE(aux)
  !
  RETURN 
  !
END SUBROUTINE


!-------------------------------------------------------------------------------
SUBROUTINE sternheimer_kernel_old(first_iter, npert, i_ref, lrdvpsi, iudvpsi, &
           thresh, dvscfins, all_conv, averlt, drhoscf, dbecsum, delta_vg)
  !-----------------------------------------------------------------------------
  !
  USE control_kcw
  USE apply_DPot_mod,        ONLY : apply_DPot_bands
  USE qpoint,                ONLY : npwq, nksq, ikks, ikqs
  USE wvfct,                 ONLY : npw, npwx, nbnd, et
  USE wavefunctions,         ONLY : evc
  USE eqv,                   ONLY : dvpsi, DPsi, evq
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE control_lr,            ONLY : lgamma, nbnd_occ
  USE buffers,               ONLY : save_buffer, get_buffer
  USE uspp_init,             ONLY : init_us_2
  USE mp,                    ONLY : mp_sum
  USE fft_base,              ONLY : dfftp, dffts
  USE klist,                 ONLY : xk, wk, ngk, igk_k
  USE uspp,                  ONLY : vkb
  USE ions_base,             ONLY : nat
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : npol
  USE units_lr,              ONLY : iuwfc, lrwfc, iudwf, lrdwf 
  USE io_global,             ONLY : stdout
  USE mp_global,             ONLY : inter_pool_comm 
  USE gvecs,                 ONLY : ngms
  USE incdrhoscf_mod,        ONLY : incdrhoscf
  !
  LOGICAL, INTENT (IN) :: first_iter
  INTEGER, INTENT (IN) :: i_ref
  INTEGER :: ik, ikk, ikq
  INTEGER, INTENT(IN) :: npert
  REAL(DP), INTENT(IN) :: thresh 
  COMPLEX(DP), POINTER, INTENT(IN) :: dvscfins(:, :, :)
  LOGICAL, INTENT(OUT) :: all_conv
  REAL(DP), INTENT(OUT) :: averlt
  COMPLEX(DP), INTENT(INOUT) :: drhoscf(dffts%nnr, nspin, npert)
  COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin, npert)
  INTEGER, INTENT(IN) :: lrdvpsi
  INTEGER, INTENT(IN) :: iudvpsi
  COMPLEX(DP) , ALLOCATABLE :: aux2(:, :)
  LOGICAL :: conv_root
  REAL(DP) :: aux_avg (2)
  REAL(DP) , ALLOCATABLE   :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian
  INTEGER :: ipert, nrec
  EXTERNAL ch_psi_all,  cg_psi
  INTEGER :: lter, ltaver, lintercall, ibnd, spin_ref
  REAL(DP) :: anorm, weight
  !
  !## DEBUG
  COMPLEX(DP) :: drhok(dffts%nnr)
  COMPLEX(DP) :: delta_vg(ngms,nspin)
  !##DEBUG 
  ! 
  all_conv = .TRUE.
  ALLOCATE (aux2(npwx*npol, nbnd))
  ALLOCATE( h_diag(npwx, nbnd) )
  ltaver = 0
  !
  linterCALL = 0
  !$acc enter data create(aux2(1:npwx*npol, 1:nbnd), dvpsi, h_diag) copyin(dpsi)
  ! 
  do ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     if (lsda) current_spin = isk (ikk)
     !
     ! read unperturbed wavefunctions psi(k) and psi(k+q)
     !
     if (nksq.gt.1) then
        if (lgamma) then
           CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
        else
           CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
           CALL get_buffer (evq, lrwfc, iuwfc, ikq)
           !$acc update device(evq)
        endif
     endif
     !
     ! compute beta functions and kinetic energy for k-point ikq
     ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
     !
     CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
     CALL g2_kin (ikq) 
     !
     ! compute preconditioning matrix h_diag used by cgsolve_all
     !
     CALL h_prec (ik, evq, h_diag)
     !
     ipert = 1
     nrec = (ipert - 1) * nksq + ik
     !
     CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, ik)
     !$acc update device(dvpsi)
     !
     ! compute the right hand side of the linear system due to
     ! the perturbation, dvscfin used as work space
     !
     IF (.NOT. first_iter ) THEN 
        ! 
        !$acc kernels present(aux2)
        aux2(:,:) = (0.D0,0.D0)
        !$acc end kernels
        CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvscfins(:, :, ipert), evc, aux2)
        !$acc kernels present(aux2, dvpsi)
        dvpsi = dvpsi + aux2
        !$acc end kernels
        !
        !
    ENDIF
     !
     ! iterative solution of the linear system (H-eS)*DPsi=dvpsi,
     ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
     !
     ! Ortogonalize dvpsi to valence states: ps = <evq|dvpsi>
     ! Apply -P_c^+. 
     ! And finally |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
     !
     CALL orthogonalize(dvpsi, evq, ikk, ikq, DPsi, npwq, .false.)
     IF ( first_iter ) THEN
        !
        !  At the first iteration DPsi is set to zero
        !
        !$acc kernels present(dpsi)
        DPsi(:, :) = (0.d0,0.d0)
        !$acc end kernels
     ELSE
        !
        ! starting value for delta_psi is read from iudwf
        !
        CALL get_buffer(DPsi, lrdwf, iudwf, nrec)
        !$acc update device(dpsi)
      ENDIF
     !
     weight = wk (ikk); conv_root = .true.
     CALL cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, DPsi, &
                       h_diag, npwx, npwq, thresh, ik, lter, conv_root, &
                       anorm, nbnd_occ(ikk), npol)
     !
     IF (.NOT. conv_root) THEN
        all_conv = .FALSE.
        WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
           &root not converged, thresh < ', es10.3)") ik, anorm
     ENDIF
     ! writes delta_psi on iunit iudwf, k=kpoint,
     !
     !$acc update self(dpsi)
     CALL save_buffer(DPsi, lrdwf, iudwf, nrec)
     !
     ltaver = ltaver + lter
     linterCALL = linterCALL + 1
     if (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4,  &
          &              " solve_linter_iu: root not converged ",e10.3)') &
          &              ik , ibnd, anorm
       
     !
     ! ...            Calculates drho, sum over k             ...
     ! ... Eventually freeze the orbital i_ref given in input ...
     !
     IF (kcw_at_KS .AND. fix_orb) THEN
       IF (spin_ref == current_spin) THEN 
          !$acc kernels present(dpsi)
          DPsi(:,i_ref) = (0.D0, 0.D0)  !! disregard the variation of the orbital i_ref 
          !$acc end kernels 
       ENDIF
     ENDIF
     !
     !$acc data copy(dbecsum(:,:,current_spin, 1) )
     CALL incdrhoscf (drhoscf(1,current_spin,1), weight, ik, &
                      dbecsum(1,1,current_spin,1), DPsi)
     !$acc end data 
    !
!!## DEBUG 
!   !! This is to check which k points are effectively equivalent (DEBUG/UNDERSTAND)
!   !! only sum over the band
!   drhok = CMPLX(0.D0,0.D0,kind=DP)
!   CALL incdrhoscf ( drhok, weight, ik, &
!                     dbecsum(1,1,current_spin,1), DPsi)
!   aux(:) = drhok(:)
!   CALL fwfft ('Rho', aux, dffts)
!   !! drhog_scf used as workspace
!   drhog_scf(:,1) = aux(dffts%nl(:))
!   !!WRITE(*,'("NICOLA", 6f22.18)') drhok(1:3) 
!   !!WRITE(*,'("NICOLA", 6f22.18)') DPsi(1:3,1) 
!   !!WRITE(*,'("NICOLA", 6f22.18)')  drhog_scf(1:3,1)
!   !!WRITE(*,'("NICOLA", 6f22.18)')  delta_vg(1:3,1)
!   WRITE(*,'("NICOLA", i5, f12.4, 3x, 2f12.8)') ik, weight, sum (CONJG(drhog_scf(:,1)) * delta_vg(:,1))*omega
!!## DEBUG
    !
  enddo ! on k-points
  !$acc exit data delete(aux2, h_diag) copyout(dvpsi, dpsi)
  !
#ifdef __MPI
   aux_avg (1) = DBLE (ltaver)
   aux_avg (2) = DBLE (lintercall)
   CALL mp_sum ( aux_avg, inter_pool_comm )
   averlt = aux_avg (1) / aux_avg (2)
#else
   averlt = DBLE (ltaver) / lintercall
#endif
   ! 
   DEALLOCATE (h_diag, aux2)
   !
END SUBROUTINE

END MODULE
