!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE rho_of_q (rhowann, ngk_all, igk_k_all, iq, rhowann_k)
  !-----------------------------------------------------------------------
  !
  !! This sunroutine compute the periodic part of the Wannier density 
  !! rho(r) = \sum_q exp[iqr]rho_q(r) 
  !! rho_q(r) = \sum_k u^*_{k,n}(r) u_{k+q,n}(r)
  !! The k+q is rewritten as p+G with p in the 1BZ, then u_{k+q,n}(r) = u_p(r)exp{-iGr}
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE fft_base,             ONLY : dffts
  USE klist,                ONLY : nkstot, xk, igk_k, ngk, nks, wk
  USE mp,                   ONLY : mp_sum
  USE control_kcw,          ONLY : evc0, iuwfc_wann, iuwfc_wann_allk, spin_component, num_wann, &
                                   mp1, mp2, mp3
  USE buffers,              ONLY : get_buffer, save_buffer
  USE wvfct,                ONLY : npwx, wg
  USE noncollin_module,     ONLY : npol
  USE control_kcw,          ONLY : map_ikq, shift_1bz, kcw_iverbosity, irr_bz
  USE cell_base,            ONLY : at
  USE mp_pools,             ONLY : inter_pool_comm
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  ! USE mp_world,             ONLY : mpime
  ! 
#ifdef DEBUG
  USE gvect,                ONLY : g, ngm
#endif
  !
  IMPLICIT NONE
  ! 
  INTEGER :: ik, ikq, npw_k, npw_kq
  !! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
  !
  INTEGER, INTENT(IN) :: iq
  !! index of the q point at hand
  !
  INTEGER :: iband, lrwfc 
  !! Band counter
  !
  COMPLEX(DP), INTENT(OUT) :: rhowann(dffts%nnr, num_wann), rhowann_k(dffts%nnr, nkstot/nspin, num_wann)
  !! The periodic part of the wannier orbital density
  !
  COMPLEX(DP) ::  evc0_kq(npwx*npol, num_wann)
  !! Auxiliary vector to store the wfc at k+q
  !
  REAL(DP) :: g_vect(3), xk_(3)
  !! G vector that shift the k+q inside the 1BZ
  !
  COMPLEX(DP) :: evc_k_g (npwx*npol), evc_k_r (dffts%nnr), phase(dffts%nnr)
  !! Auxiliary wfc in reciprocal and real space, the phase associated to the hift k+q-> k'
  !
  COMPLEX(DP) :: evc_kq_g (npwx*npol), evc_kq_r (dffts%nnr)
  !! Auxiliary wfc in reciprocal and real space
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !! The global index of k-points
  !
  INTEGER :: global_ik, ik_eff
  !
  INTEGER, INTENT(IN) :: &
       igk_k_all(npwx,nkstot),&    ! index of G corresponding to a given index of k+G
       ngk_all(nkstot)             ! number of plane waves for each k point
  !
  INTEGER :: isym
  REAL(DP) :: weight
  !
#ifdef DEBUG
  INTEGER :: ig_save, ig, ip, ir
  COMPLEX (DP ) :: pippo
  REAL(DP) :: pippo_real, xq(3)
#endif
  !  
  CALL start_clock ( 'rho_of_q' )
  DO ik = 1, nks
    ! CHECK: Need to understand/think more about pool parallelization
    ! what happen if k+q is outside the pool??
    ! Some problem in EXX (have a look to /PW/src/exx.f90 exx_grid_init)
    ! Solved for now with WFs Broadcast (See bcast_wfc.f90) 
    !
    IF ( lsda ) current_spin = isk(ik)
    IF ( lsda .AND. current_spin /= spin_component) CYCLE
    !
    weight = 1.D0/(mp1*mp2*mp3)
    !
    xk_ = xk(:,ik)
    CALL cryst_to_cart(1, xk_, at, -1)
    !
#ifdef DEBUG
     WRITE(*,'(10x, "ik = ", i5, 3x, "xk = ", 3f12.6, " [Cryst]")') ik, (xk_(ip), ip=1,3) 
#endif
    !
    global_ik = global_kpoint_index (nkstot,ik)
    global_ik = global_ik - (spin_component -1)*nkstot/nspin
    lrwfc = num_wann*npwx
    ik_eff = ik-(spin_component-1)*nkstot/nspin
    CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
    !! ... Retrive the ks function (in the Wannier Gauge)
    !
    !IF (irr_bz) THEN 
    !   CALL find_index_1bz_smart_ibz(global_ik, iq, g_vect, ikq, isym) 
    !ELSE
       CALL find_index_1bz_smart(global_ik, iq, g_vect, ikq) 
    !ENDIF 
    ! see also other implementations (more robust?)of the search inside compute_map_ikq_single.f90 
    ! if the search fails
    !
    !ikq = map_ikq(global_ik)
    !g_vect(:) = shift_1bz(:,global_ik)
    !! ... The index ikq in the 1BZ corresponding at the k+q, and the vector G_bar defining the shift 
    !! ... xk(:,ikq)+G_bar = xk(:,k+q)
    !! ... see compute_map_ikq_single.f90
    !
#ifdef DEBUG
    !! WORKS ONLY on 1 proc!! 
    ig_save = 0
    do ig =1 ,ngm 
       ! in general for k/=gamma the index corresponding to a given g vector is different from the one at k=gamma.
       ! here I scan the gvectors k+g (given by igk(i,ikq) ) untill I found the one that matches with the cohoosen 
       ! one g_vect(:), and I store the index to be used after. This is the index ii in the evc_k(ii) corresponding 
       ! to the g vector g_vect  
       ! igk_k(i,ik) is the index in the original list of g  vectors (the one at GAMMA) of the i-th component 
       ! in evc_k(i) See PW/src/pwcomm, klist module.
       pippo_real = sqrt( sum( (g_vect(:)-g(:,igk_k(ig,ikq)))**2 ) )
       ! write(*,'(2i5, f12.6, 3x, 3f12.6, 3x, 3f12.6)') ikq, ig, pippo_real, g_vect(:), g(:,igk_k(ig,ikq))
       if ( abs(pippo_real) .lt. 1e-5) then 
         ig_save = ig
         xq = g_vect 
         CALL cryst_to_cart(1, xq, at, -1)
         WRITE(*,'(10x, "DEBUG: MATCH FOUND for g=",3f12.6, " [Cryst]", 3x, "index_k+g = ", i3, " -->", i3, /)') (xq(ir),ir=1,3), ig_save, igk_k(ig_save,ikq)
         exit
       endif
    enddo
    IF( ig_save == 0 ) CALL errore('do_koopmans', 'G vector shift NOT FOUND in the list of ngm g_vectors', 1)
#endif
    !
    phase(:) = 0.D0
    CALL calculate_phase(g_vect, phase) 
    !! ... Calculate the phase associated to the k+q-> ikq map: exp[ -i(G_bar * r) ]
    !
    lrwfc = num_wann * npwx 
    CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann_allk, ikq )
    !! ... Retrive the ks function (in the Wannier Gauge) at the k-point in the 1BZ
    !! ... (possibly in the Irreducible wedge) connected to k+q 
    !
    DO iband = 1, num_wann
       !
       !IF (irr_bz) weight = wk(ik)*nspin
!       WRITE(*,'("ik, ibndn, weight=", 2I5, 3F12.6)') ik, iband, weight 
       !
       npw_k = ngk(ik)
       evc_k_g(:) =  evc0(:,iband)
       evc_k_r(:) = ZERO
       CALL invfft_wave (npw_k, igk_k (1,ik), evc_k_g , evc_k_r )
       !! ... The wfc in R-space at k
       !
       npw_kq = ngk_all(ikq)
       evc_kq_g = evc0_kq(:,iband)
       evc_kq_r = ZERO
       CALL invfft_wave (npw_kq, igk_k_all (1,ikq), evc_kq_g , evc_kq_r )
       !
       !IF (irr_bz) THEN 
       !   CALL rotate_evc(isym, phase, evc_kq_r)
       !ELSE 
          evc_kq_r(:) = evc_kq_r(:)*phase(:)
       !ENDIF 
       !! ... Rotate the KS state in the IBZ with the symmetry operation to get the KS state
       !! ... at k+q 
       !
       ! ... The wfc in R-space at S(k')+G_bar <-- k+q where S(k')-G_bar = (k+q)
       ! ... evc_k+q(r) = sum_G exp[iG r] c_(k+q+G) = sum_G exp[iG r] c_S(k')+G_bar+G 
       !            = exp[-iG_bar r] sum_G' exp[iG'r] c_S(k')+G' = exp[-iG_bar r] *evc_S(k')(r) = exp[-iG_bar r] *evc_k'(S^-1(r))
       !
       !rhowann(:,iband) = rhowann(:,iband) + conjg(evc_k_r(:))*evc_kq_r(:)*phase(:)/weight !*wg(iband,ik)
       rhowann(:,iband) = rhowann(:,iband) + conjg(evc_k_r(:))*evc_kq_r(:)*weight !*wg(iband,ik)
       rhowann_k(:,ik,iband) = conjg(evc_k_r(:))*evc_kq_r(:)*weight !*wg(iband,ik)
       ! ... The periodic part of the wannier-orbital density in real space
       ! ... rho_q(r) = sum_k [ evc_k,v(r)* evc_k+q,v(r)] = sum_k [ evc_k,v(r)* evc_k',v(r) exp[-iG_bar r]]
       ! 
#ifdef DEBUG
       ! This is to check that the phase factor is correct. 
       ! pippo is the Fourier component at g_vect of evc_kq(r) computed explicitely (without calling FFT)  
       ! It has to match with the g_vect component of evc_kq(g) whose index ig_save was found above
       ! If phase was correctly computed pippo-evc_kq(ig_save) HAS TO BE IDENTICALLY ZERO
       pippo = (0.D0,0.d0)
       DO ir = 1, dffts%nnr
         pippo = pippo + evc_kq_r(ir)*phase(ir)
       ENDDO
       pippo = pippo/dffts%nnr
       pippo_real = REAL(pippo-evc_kq_g(ig_save))**2+AIMAG(pippo-evc_kq_g(ig_save))**2
       IF ( pippo_real .gt. 1e-6) THEN 
          WRITE(*, '("ikq =" i5, 3x, "nbnd =" i5, 3x,  "CHECK =" 9f15.8, 2i5 )') & 
          ikq, iband, pippo, evc_kq_g(ig_save), pippo-evc_kq_g(ig_save), g(:,igk_k(ig_save,ikq)), igk_k(ig_save,ikq), ig_save
          CALL errore('do_koopmans','Something WRONG with the pahse factor',1)
       ENDIF
#endif 
      !
    ENDDO ! bands
    ! 
  ENDDO ! kpoints
  !
  CALL mp_sum( rhowann, inter_pool_comm )
  ! ... Sum up the results over k point in different pools
  !
  CALL stop_clock ('rho_of_q')
  !
END subroutine
