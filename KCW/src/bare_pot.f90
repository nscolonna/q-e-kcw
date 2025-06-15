! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!#define DEBUG
!!JA Assume all array are on GPU
!-----------------------------------------------------------------------
SUBROUTINE bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
  !---------------------------------------------------------------------
  !
  !! This routine compute the Hxc potential due to the (peridoc part of the) wannier
  !! charge density. V^{0n}_{Hxc}(r) = \int dr' f_Hxc(r,r') w^{0n}(r')
  !  
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvecs,                ONLY : ngms
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : tpiba2, omega
  USE control_kcw,          ONLY : spin_component, kcw_iverbosity, x_q, nrho
  USE gvect,                ONLY : g
  USE qpoint,               ONLY : xq
  USE constants,            ONLY : e2, fpi
  USE eqv,                  ONLY : dmuxc
  USE control_lr,           ONLY : lrpa
  USE martyna_tuckerman,    ONLY : wg_corr_h, do_comp_mt
  USE io_global,            ONLY : stdout
  !USE exx_base,             ONLY : g2_convolution
  USE coulomb,             ONLY : g2_convolution
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_lsda, nspin_gga, nspin_mag, npol
  ! GC LR suff
  USE xc_lib,            ONLY : xclib_dft_is
  USE scf,               ONLY : rho, rho_core
  USE uspp,              ONLY : nlcc_any
  USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE dv_of_drho_lr,     ONLY : dv_of_drho_xc


  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT (OUT) :: rhog(ngms,nrho)
  ! ... periodic part of wannier density in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg(ngms,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vg_(ngms,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in G-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (OUT) :: vh_rhog(ngms)
  ! ... Hartree perturbing potential [f_hxc(r,r') x wann(r')] in G-space
  ! 
  COMPLEX(DP), INTENT (OUT) :: delta_vr(dffts%nnr,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space
  !
  COMPLEX(DP), INTENT (OUT) :: delta_vr_(dffts%nnr,nspin_mag)
  ! ... perturbing potential [f_hxc(r,r') x wann(r')] in r-space without g=0 contribution
  !
  COMPLEX(DP), INTENT (IN)  :: rhor(dffts%nnr,nrho)
  ! ... periodic part of wannier density in r-space
  !
  INTEGER, INTENT (IN)      :: iq
  ! ... q-point index
  !
  !!COMPLEX(DP)               :: aux(dffts%nnr), aux_(dffts%nnr)
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:)        :: aux, aux_
  ! ... auxiliary vectors 
  !
  COMPLEX(DP), ALLOCATABLE  :: vaux(:)
  ! ... auxiliary vector
  !
  INTEGER                   :: ig, is, ir, iss, ip
  ! ... counters 
  !
  REAL(DP)                  :: qg2, eh_corr
  ! ... |q+G|^2, g=0 correction to the hartree energy 
  !
  LOGICAL                   :: lgamma
  !
  REAL(DP)                  :: xkq(3), xk(3)
  ! ... coordinate of k and k+q 
  !
  !!REAL(DP)                  :: fac(ngms)
  REAL(DP), ALLOCATABLE   :: fac(:)
  ! ... Coulomb kernel possibly with the special treatment of the q+g+0 component 
  !
  !!COMPLEX(DP)               :: vh_rhog_g0eq0(ngms)
  COMPLEX(DP), ALLOCATABLE   :: vh_rhog_g0eq0(:)
  ! ... The hartree potential with th q+g=0 component set to zero
  !
  COMPLEX(DP), ALLOCATABLE :: rhor_(:,:)
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nls_d(:)
#else
  INTEGER, ALLOCATABLE :: nls_d(:)
#endif


#if defined(__CUDA)
  nls_d  => dffts%nl_d
#else
  ALLOCATE( nls_d(dffts%ngm) )
  nls_d  = dffts%nl
#endif
  !
  ALLOCATE(aux(dffts%nnr))
  !$acc enter data create(aux)
  !$acc data present(rhor, rhog, vh_rhog, delta_vr, delta_vg, delta_vr_, delta_vg_)
  !! The periodic part of the orbital density in g space  
  DO ip=1,nrho !<---CONSIDER TO SUBSTITUTE WITH nspin_mag
      !$acc kernels present(aux)
      aux(:) = rhor(:,ip)/omega
      !$acc end kernels
      !$acc host_data use_device(aux)
      CALL fwfft ('Rho', aux, dffts) 
      !$acc end host_data
      !$acc kernels present(aux) deviceptr(nls_d)
      rhog(:,ip) = aux(nls_d(:))
      !$acc end kernels
  END DO
  !
 !!JA aux      = ZERO
  !
  ! .. First the xc contribution
  !
  IF (.NOT. lrpa) THEN 
    ALLOCATE ( rhor_(dffts%nnr,nspin_mag) )
!GPU    !$acc enter data create(rhor_)
!GPU    !$acc kernels present(rhor_)
    rhor_ = CMPLX(0.D0,0.D0,kind=DP)
!GPU    !$acc end kernels
    IF (nspin_mag == 4) THEN
!GPU      !$acc kernels present(rhor_)
      rhor_=rhor/omega
!GPU      !$acc end kernels
    ELSE
!GPU      !$acc kernels present(rhor_)
      rhor_(:,spin_component) = rhor(:,1)/omega
!GPU      !$acc end kernels
    ENDIF
!GPU    !$acc enter data copyin(dmuxc) 
!GPU    !$acc kernels 
      delta_vr = (0.0_dp, 0.0_dp)
!GPU    !$acc end kernels  
    CALL dv_of_drho_xc(delta_vr, rhor_) !!JA ON GPU
    !$acc update device(delta_vr)
!GPU    !$acc exit data delete(rhor_, dmuxc)
    DEALLOCATE (rhor_)
    !$acc kernels 
    delta_vr_ = delta_vr 
    !$acc end kernels
  ELSE
    !$acc kernels 
    delta_vr = ZERO
    delta_vr_ = ZERO
    !$acc end kernels
  ENDIF 
  !
  !
  xk(:) = 0.D0
  xkq(:) = -x_q(:,iq)
  ! ... auxiliary coordinate to pass to g2_convoltion 
  ! 
  ! ... The Hartree contribution first 
  !
  ALLOCATE(fac(ngms))
  CALL g2_convolution(ngms, g, xk, xkq, fac)
  ! ... the hartree kernel (eventually within the 
  ! ... Gygi-Baldereschi sheme, see setup_coulomb) 
  !
  ALLOCATE(vh_rhog_g0eq0(ngms))
  !$acc enter data create(vh_rhog_g0eq0) copyin(fac)
  !$acc parallel loop private(qg2) present(fac, vh_rhog_g0eq0) present_or_copyin(g, x_q(:,iq))
  DO ig = 1, ngms
    !
   !! qg2 = SUM ( (g(:,ig)+x_q(:,iq))**2 )
    qg2 =  (g(1,ig)+x_q(1,iq))**2  +  (g(2,ig)+x_q(2,iq))**2  +  (g(3,ig)+x_q(3,iq))**2  
    !
    IF (qg2 .lt. 1.0e-8_dp) THEN 
        vh_rhog_g0eq0(ig) = (0.0_dp, 0.0_dp)
    ELSE
        vh_rhog_g0eq0(ig) =  e2 * fpi * rhog(ig,1) / (tpiba2 * qg2)
    END IF
    ! ... set to zero the q+g=0 component
    !
    vh_rhog(ig) =  rhog(ig,1) * cmplx(fac(ig), 0.d0, KIND=dp)
    ! ... the hartree potential possibly with the special treatment of the q+g=0 component  
    !
  ENDDO
  !$acc exit data delete(fac)
  DEALLOCATE(fac)
  !
  ! ... eventually add MT corrections
  !
  lgamma = (xq(1)==0.D0.AND.xq(2)==0.D0.AND.xq(3)==0.D0)
  IF (do_comp_mt .AND. lgamma) then
     ! ... this make sense only for a GAMMA only calculation 
     ! ... (do_compt_mt = .false if more than one k-point (see kcw_readin.f90) 
     !
     IF (kcw_iverbosity .gt. 1 ) WRITE(stdout,'(5x, " ADDING Martyna-Tuckerman correction" ,/)')
     !
     ALLOCATE( vaux( ngms ) )
     !$acc update self(rhog)
     CALL wg_corr_h (omega, ngms, rhog, vaux, eh_corr)
     !$acc enter data copyin(vaux)
     !$acc kernels
     vh_rhog(1:ngms) = vh_rhog(1:ngms) +  vaux(1:ngms)
     vh_rhog_g0eq0(1:ngms) = vh_rhog_g0eq0(1:ngms) +  vaux(1:ngms)
     !$acc end kernels
     !$acc exit data delete(vaux)
     DEALLOCATE( vaux )
  ENDIF
  !
  ! ... Go to r-space 
  !
  ALLOCATE(aux_(dffts%nnr))
  !$acc enter data create(aux_)
  !$acc kernels present(aux, aux_)
  aux=(0.d0,0.d0)
  aux_=(0.d0,0.d0)
  !$acc end kernels
  !$acc kernels present(aux, aux_) deviceptr(nls_d)
  aux(nls_d(:))  = vh_rhog(:)                    ! G-space components of the Hartree potential
  aux_(nls_d(:)) = vh_rhog_g0eq0(:)
  !$acc end kernels
  !$acc exit data delete(vh_rhog_g0eq0)
  DEALLOCATE(vh_rhog_g0eq0)
  !$acc host_data use_device(aux)
  CALL invfft ('Rho', aux, dffts)
  !$acc end host_data
  !$acc host_data use_device(aux_)
  CALL invfft ('Rho', aux_, dffts)
  !$acc end host_data
  !
  IF (nspin_mag==4 .and. domag) THEN    ! Perturbing potential due to Hartree
    !$acc kernels present(aux, aux_)
    delta_vr(:,1) = delta_vr(:,1)   + aux(:)
    delta_vr_(:,1) = delta_vr_(:,1) + aux_(:)
    !$acc end kernels
  ELSEIF (nspin_mag==4 .and. .NOT. domag) THEN
    !$acc kernels present(aux, aux_)
    delta_vr(:,1) = delta_vr(:,1) + aux(:)
    delta_vr_(:,1) = delta_vr_(:,1) + aux_(:)
    !$acc end kernels
  ELSE
    !$acc kernels present(aux, aux_)
    DO is = 1, nspin_mag
      delta_vr(:,is)  = delta_vr(:,is)   + aux(:)
      delta_vr_(:,is) = delta_vr_(:,is) + aux_(:)
    END DO
    !$acc end kernels
  END IF
  !
  DO is = 1, nspin_mag
    !
    !$acc kernels present(aux, aux_)
    aux(:) = delta_vr(:,is)
    aux_(:) = delta_vr_(:,is) 
    !$acc end kernels
    !
    !$acc host_data use_device(aux)
    CALL fwfft ('Rho', aux, dffts)
    !$acc end host_data
    !$acc host_data use_device(aux_)
    CALL fwfft ('Rho', aux_, dffts)
    !$acc end host_data
    !
    !$acc kernels present(aux, aux_) deviceptr(nls_d)
    delta_vg(:,is)  = aux(nls_d(:))
    delta_vg_(:,is) = aux_(nls_d(:))
    !$acc end kernels
    !
  ENDDO
  !$acc exit data delete(aux, aux_)
  !$acc end data
  DEALLOCATE(aux, aux_)
  !
  !
END subroutine bare_pot
