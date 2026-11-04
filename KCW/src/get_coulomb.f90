
SUBROUTINE write_coulomb()
USE kinds,                ONLY : DP
USE control_kcw,          ONLY : Vcoulomb, Wcoulomb
USE control_kcw,          ONLY : spin_component, get_coulomb, irvect_shifted
USE control_kcw,          ONLY : num_wann, nqstot, num_R
USE noncollin_module,     ONLY : nspin_mag
!
IMPLICIT NONE
!
INTEGER              ::    iun_coulomb
CHARACTER(len=1024)  ::    filename  
INTEGER              ::    iwann, jwann, ir, is
INTEGER              ::    spin_index

!
!
  ! Vcoulomb/Wcoulomb(is,...) are stored by spin_index (coulomb_me): is = 1 same
  ! spin as spin_component, is = 2 opposite spin. Write is itself as the label:
  ! spin_index(is) swapped the two labels when spin_component = 2.
  !
  filename = 'barecoulomb.txt'
  iun_coulomb = 237 
  OPEN (iun_coulomb, file = filename)

  DO ir = 1, num_R
    WRITE(iun_coulomb, '(3I7)') irvect_shifted(:, ir) 
    DO iwann=1, num_wann
      DO jwann=1, num_wann
        DO is = 1, nspin_mag !one for spin component, other for non spin component
          WRITE(iun_coulomb, '(3I5, 3X, 2ES24.16)') iwann, jwann, is, &
          real(Vcoulomb(is, ir, jwann, iwann)), aimag(Vcoulomb(is, ir, jwann, iwann)) 
        END DO!is
      END DO!jwann
    END DO!iwann
  END DO!ir
  CLOSE(iun_coulomb)
  !
  !
  iun_coulomb = 238 
  filename = 'screencoulomb.txt'
  OPEN (iun_coulomb, file = filename)
  !
  DO ir = 1, num_R
    WRITE(iun_coulomb, '(3I7)') irvect_shifted(:, ir) 
    DO iwann=1, num_wann
      DO jwann=1, num_wann
        DO is = 1, nspin_mag !one for spin component, other for non spin component
          WRITE(iun_coulomb, '(3I5, 3X, 2ES24.16)') iwann, jwann, is, &
               real(Wcoulomb(is, ir, jwann, iwann)), aimag(Wcoulomb(is, ir, jwann, iwann)) 
        END DO!is
      END DO!jwann
    END DO!iwann
  END DO!ir
  CLOSE(iun_coulomb)
!
END SUBROUTINE write_coulomb

SUBROUTINE coulomb_me( iwann, iq, drhog_scf, rhowann, weight_q)

!this function calculates the Coulomb interaction matrix element:
! <rho_{R,iwann}|V_{Hxc}|rho_{0,jwann}>
!
USE kinds,                ONLY : DP
USE control_kcw,          ONLY : nqstot, tmp_dir_save, num_wann, nrho
USE control_kcw,          ONLY : iurho_wann, x_q, irvect_shifted
USE control_kcw,          ONLY : Vcoulomb, Wcoulomb
USE control_kcw,          ONLY : spin_component, get_coulomb, num_R
USE buffers,              ONLY : get_buffer 
USE io_files,             ONLY : tmp_dir
USE fft_base,             ONLY : dffts
USE klist,                ONLY : nkstot
USE gvecs,                ONLY : ngms
USE noncollin_module,     ONLY : nspin_mag
USE constants,            ONLY : tpi
USE lsda_mod,             ONLY : nspin
USE cell_base,            ONLY : omega, at
USE mp_bands,             ONLY : intra_bgrp_comm
USE mp,                   ONLY : mp_sum
!
IMPLICIT NONE
COMPLEX(DP)                  :: IMAG = (0.D0,1.D0)
INTEGER, INTENT(IN)          :: iwann
! (fixed) index of wannier function
INTEGER, INTENT(IN)          :: iq
! (fixed) index of q point
! quantity V_xc|rho_i> in G space
COMPLEX(DP), INTENT(IN)      :: drhog_scf(ngms,nspin_mag)
!LR wannier density for fixed q and iwann
COMPLEX(DP), INTENT(IN)      :: rhowann(dffts%nnr, num_wann, nrho)
INTEGER                      :: jwann
! index of wannier functions to loop over
REAL(DP), INTENT(IN)         :: weight_q 
INTEGER                      :: ip
! index of nrho to loop over
INTEGER                      :: ir
!
!quantities needed for computing bare_pot for jwann
!
COMPLEX(DP), ALLOCATABLE     :: rhog(:,:), rhog_iwann(:,:)
! ... periodic part of wannier density in G-space
COMPLEX(DP), ALLOCATABLE     :: delta_vg(:,:)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space
COMPLEX(DP), ALLOCATABLE     :: delta_vg_(:,:)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space without g=0 contribution
COMPLEX(DP), ALLOCATABLE     :: vh_rhog(:)
! ... Hartree perturbing potential [f_hxc(r,r') x wann(r')] in G-space
COMPLEX(DP), ALLOCATABLE     :: delta_vr(:,:)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space
COMPLEX(DP), ALLOCATABLE     :: delta_vr_(:,:)
! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space without g=0 contribution
COMPLEX(DP), ALLOCATABLE     :: rhor(:,:)
! ... periodic part of wannier density in r-space
INTEGER                      :: lrrho
INTEGER                      :: is, is_, is1, js, spin_save
INTEGER                      :: spin_index
COMPLEX(DP), ALLOCATABLE     :: delta_vg_js(:,:,:), delta_vg_js_(:,:,:)
! ... perturbing potential of rho_jwann placed in spin channel js (last index),
!     with and without the q+G=0 component
!
! we are already inside loops over iq and iwann
!
COMPLEX(DP)                  :: pi_q_unrelax(2), pi_q_relax(2), pi_q_unrelax_(2)
REAL(DP)                     ::x_q_cryst(3)
!
!
ALLOCATE(rhog(ngms,nrho))
ALLOCATE(rhog_iwann(ngms,nrho))
ALLOCATE(delta_vg(ngms,nspin_mag))
ALLOCATE(delta_vg_(ngms,nspin_mag))
ALLOCATE(vh_rhog(ngms))
ALLOCATE(delta_vr(dffts%nnr,nspin_mag))
ALLOCATE(delta_vr_(dffts%nnr,nspin_mag))
ALLOCATE(rhor(dffts%nnr,nrho))
ALLOCATE(delta_vg_js(ngms,nspin_mag,nspin_mag))
ALLOCATE(delta_vg_js_(ngms,nspin_mag,nspin_mag))
!
IF (nrho==4) THEN
    CALL errore('output_coulomb', 'non-collinear not implemented &
      for coulomb matrix elements.', 1)
END IF
!
ip=1
rhor(:,ip) = rhowann(:,iwann,ip)
CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_) 
rhog_iwann= rhog
!
DO jwann = 1, num_wann
  !here rhowann is already filled with the wannier density in real space 
  !for fixed iq
  !
  ! get rhog for jwann by fourier transforming rhor 
  !WARNING! WON'T WORK FOR NON_COLLINEAR
    DO ir = 1, num_R
      ip=1
      rhor(:,ip) = rhowann(:,jwann,ip)
      ! rho_q(r-R) = rho_q(r) e^{-iqR}
      x_q_cryst(:) = x_q(:,iq)
      CALL cryst_to_cart(1,x_q_cryst,at,-1)
      rhor(:, ip) = rhor(:, ip) * EXP( -IMAG*tpi*DOT_PRODUCT(x_q_cryst(:),irvect_shifted(:,ir)) )
      !
      ! Perturbing potential of rho_m(R) (m = jwann) placed in each spin channel js:
      !   delta_vg_js(:,s1,js) = [ v + f_xc^{s1,js} ] rho_m(R)       (G-space)
      ! bare_pot puts the density in the channel spin_component: switch it
      ! temporarily (it is restored before spin_index is used below).
      !
      spin_save = spin_component
      DO js = 1, nspin_mag
        spin_component = js
        CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg_js(:,:,js), iq, &
                        delta_vr_, delta_vg_js_(:,:,js) )
      END DO
      spin_component = spin_save
      !
      ! for now we only work with the density, i.e. ip = 1.
      ! Component spin_index(js) is the element between rho_n (n = iwann) in the
      ! channel spin_component (the one of the LR perturbation, drhog_scf) and
      ! rho_m(R) in the channel js:
      !   bare:      <rho_n,sc | v + f_xc^{sc,js} | rho_m,js(R)>
      !   screened:  bare (without q+G=0) + sum_s1 <drho_s1[rho_n,sc] | v + f_xc^{s1,js} | rho_m,js(R)>
      ! Before, the induced (screened) part used rho_m in the channel spin_component
      ! for both components, so component 2 of W was V_opposite + (W - V)_same instead
      ! of W_opposite. The bare part is unchanged, since f_xc^{sc,js} = f_xc^{js,sc}.
      !
      pi_q_unrelax (:) = (0.D0, 0.D0)
      pi_q_unrelax_(:) = (0.D0, 0.D0)
      pi_q_relax   (:) = (0.D0, 0.D0)
      DO js = 1, nspin_mag
          pi_q_unrelax (spin_index(js)) = weight_q * omega * &
                                          SUM( CONJG(delta_vg_js (:,spin_save,js)) * rhog_iwann(:,1) )
          pi_q_unrelax_(spin_index(js)) = weight_q * omega * &
                                          SUM( CONJG(delta_vg_js_(:,spin_save,js)) * rhog_iwann(:,1) )
          !
          ! screened potential
          ! < jwann, R, js | W | 0, iwann, spin_component >
          DO is1 = 1, nspin_mag
            pi_q_relax(spin_index(js)) = pi_q_relax(spin_index(js)) + &
                 SUM( CONJG(drhog_scf(:,is1)) * delta_vg_js(:,is1,js) ) * weight_q * omega
          END DO
      END DO
      CALL mp_sum (pi_q_unrelax,  intra_bgrp_comm)
      CALL mp_sum (pi_q_unrelax_, intra_bgrp_comm)
      Vcoulomb(:, ir, jwann, iwann) = Vcoulomb(:, ir, jwann, iwann) + pi_q_unrelax(:)
      !
      CALL mp_sum (pi_q_relax, intra_bgrp_comm)
      pi_q_relax(:) = pi_q_relax(:) + pi_q_unrelax_(:)
      Wcoulomb(:, ir, jwann, iwann) = Wcoulomb(:, ir, jwann, iwann) + pi_q_relax(:)
      IF (jwann .eq. iwann .and. SUM(ABS(iRvect_shifted(:,ir))) .lt. 1.D-06 )&
      WRITE(*,*) "(TO COMPARE with output) iq = ", iq, "rvect=", irvect_shifted(:, ir), "pi_q_unrelax=", pi_q_unrelax, &
       "pi_q_relax = ", pi_q_relax
  
    END DO !ir
END DO !jwann
!
DEALLOCATE(rhog)
DEALLOCATE(rhog_iwann)
DEALLOCATE(delta_vg)
DEALLOCATE(delta_vg_)
DEALLOCATE(vh_rhog)
DEALLOCATE(delta_vr)
DEALLOCATE(delta_vr_)
DEALLOCATE(rhor)
DEALLOCATE(delta_vg_js, delta_vg_js_)


END SUBROUTINE



FUNCTION spin_index(is) 
  !! spin index for up-up / down-down -> 1 
  !! spin index for up-down /down-up  -> 2
  USE control_kcw,     ONLY : spin_component
  IMPLICIT NONE
  INTEGER    :: is 
  INTEGER    :: spin_index 

  IF( is .eq. spin_component ) THEN 
    spin_index = 1
  ELSE 
    spin_index = 2 
  END IF
END FUNCTION
