!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthogonalize_omega(dvpsi, evq, ikk, ikq, dpsi, npwq, w)
!------------------------------------------------------------------------
  !
  ! This routine orthogonalizes dvpsi to the valence states: ps = <evq|dvpsi>
  ! It should be quite general. It works for metals and insulators, with
  ! NC as well as with US PP, both SR or FR.
  ! Note that on output it changes sign. So it applies -P^+_c.
  !
  ! NB: IN/OUT is dvpsi ; dpsi is used as work_space
  !
USE kinds, ONLY : DP
USE klist, ONLY : lgauss, degauss, ngauss
USE noncollin_module, ONLY : noncolin, npol
USE wvfct, ONLY : npwx, nbnd, et
USE ener, ONLY : ef
USE control_lr,  ONLY : alpha_pv, nbnd_occ
USE becmod,      ONLY : bec_type, becp, calbec
USE uspp,        ONLY : vkb, okvan
USE mp_bands,    ONLY : intra_bgrp_comm
USE mp,          ONLY : mp_sum
USE control_flags, ONLY : gamma_only, offload_type
USE gvect,       ONLY : gstart
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: ikk, ikq   ! the index of the k and k+q points
INTEGER, INTENT(IN) :: npwq       ! the number of plane waves for q
COMPLEX(DP), INTENT(IN) :: evq(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,nbnd)
COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbnd) ! work space allocated by
                                                   ! the calling routine
COMPLEX(DP), INTENT(IN) :: w

COMPLEX(DP), ALLOCATABLE :: ps(:,:)
REAL(DP), ALLOCATABLE :: ps_r(:,:)
COMPLEX(DP), ALLOCATABLE :: wwg(:)

INTEGER :: ibnd, jbnd, nbnd_eff
REAL(DP) :: wg1, w0g, wgp, deltae, theta
REAL(DP), EXTERNAL :: w0gauss, wgauss
! functions computing the delta and theta function

CALL start_clock ('ortho')
IF (gamma_only) THEN
   ALLOCATE(ps_r(nbnd,nbnd))
   !$acc enter data create(ps_r(1:nbnd,1:nbnd))
   !$acc kernels present(ps_r)
   ps_r = 0.0_DP
   !$acc end kernels
ENDIF

ALLOCATE(ps(nbnd,nbnd))
ALLOCATE(wwg(nbnd))

!$acc data present_or_copyin(evq) present_or_copy(dvpsi,dpsi) create(ps(1:nbnd,1:nbnd))

!$acc kernels present(ps)
ps = (0.0_DP, 0.0_DP)
!$acc end kernels

!
if (lgauss) then
   !
   IF (gamma_only) CALL errore ('orthogonalize', "degauss with gamma &
        & point algorithms",1)
   !
   !  metallic case
   !
   !$acc host_data use_device(evq, dvpsi, ps)
   IF (noncolin) THEN
      CALL myzgemm( 'C', 'N', nbnd, nbnd_occ (ikk), npwx*npol, (1.d0,0.d0), &
                 evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
   ELSE
      CALL myzgemm( 'C', 'N', nbnd, nbnd_occ (ikk), npwq, (1.d0,0.d0), &
                 evq, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
   END IF
   !$acc end host_data
   !

   DO ibnd = 1, nbnd_occ (ikk)
      wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
      w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
      DO jbnd = 1, nbnd
         wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
         deltae = et (jbnd, ikq) - et (ibnd, ikk)
         theta = wgauss (deltae / degauss, 0)
         wwg(jbnd) = CMPLX(wg1 * (1.d0 - theta) + wgp * theta,0.0_DP)
         IF (jbnd <= nbnd_occ (ikq) ) THEN
            wwg(jbnd) = wwg(jbnd) + alpha_pv * theta * (wgp - wg1) / (deltae - w)
         ENDIF
      ENDDO
      !
      !$acc parallel loop copyin(wwg) present(ps)
      DO jbnd = 1, nbnd
         ps(jbnd,ibnd) = wwg(jbnd) * ps(jbnd,ibnd)
      ENDDO
      !$acc end parallel loop
      !
      IF (noncolin) THEN
         !$acc kernels present(dvpsi)
         dvpsi(1:npwx*npol, ibnd) = wg1 * dvpsi(1:npwx*npol, ibnd)
         !$acc end kernels
      ELSE
         !$acc kernels present(dvpsi)
         dvpsi(1:npwq, ibnd) = wg1 * dvpsi(1:npwq, ibnd)
         !$acc end kernels
      END IF
   END DO
   nbnd_eff=nbnd
ELSE
   !
   !  insulators
   !
   !$acc host_data use_device(evq, dvpsi, ps)
   IF (noncolin) THEN
      CALL myzgemm( 'C', 'N',nbnd_occ(ikq), nbnd_occ(ikk), npwx*npol, &
             (1.d0,0.d0), evq, npwx*npol, dvpsi, npwx*npol, &
             (0.d0,0.d0), ps, nbnd )
   ELSEIF (gamma_only) THEN
        !$acc host_data use_device(ps_r)
            CALL mydgemm( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), 2*npwq, &
             2.0_DP, evq, 2*npwx, dvpsi, 2*npwx, &
             0.0_DP, ps_r, nbnd )
            IF (gstart == 2 ) THEN
               CALL mydger( nbnd_occ(ikq), nbnd_occ (ikk), -1.0_DP, evq, &
                    & 2*npwq, dvpsi, 2*npwx, ps_r, nbnd )
            ENDIF
        !$acc end host_data
   ELSE
      CALL myzgemm( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), npwq, &
             (1.d0,0.d0), evq, npwx, dvpsi, npwx, &
             (0.d0,0.d0), ps, nbnd )
   END IF
   !$acc end host_data
   nbnd_eff=nbnd_occ(ikk)
END IF
IF (gamma_only) THEN
   !$acc host_data use_device(ps_r)
   call mp_sum(ps_r(:,:),intra_bgrp_comm)
   !$acc end host_data
ELSE
   !$acc host_data use_device(ps)
   call mp_sum(ps(:,1:nbnd_eff),intra_bgrp_comm)
   !$acc end host_data
ENDIF
!
! dpsi is used as work space to store S|evc>
!
IF (okvan) CALL calbec ( offload_type, npwq, vkb, evq, becp, nbnd_eff)
CALL s_psi (npwx, npwq, nbnd_eff, evq, dpsi)
!
! |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
!
IF (gamma_only) THEN
   !$acc kernels present(ps)
   ps = CMPLX (ps_r,0.0_DP, KIND=DP)
   !$acc end kernels
ENDIF
!
!$acc host_data use_device(dpsi, ps, dvpsi)
if (lgauss) then
   !
   !  metallic case
   !
   IF (noncolin) THEN
      CALL myzgemm( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
                (1.d0,0.d0), dpsi, npwx*npol, ps, nbnd, (-1.0d0,0.d0), &
                dvpsi, npwx*npol )
   ELSE
      CALL myzgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd, &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
              dvpsi, npwx )
   END IF
ELSE
   !
   !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
   !
   IF (noncolin) THEN
      CALL myzgemm( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd_occ(ikk), &
                (1.d0,0.d0),dpsi,npwx*npol,ps,nbnd,(-1.0d0,0.d0), &
                dvpsi, npwx*npol )
   ELSEIF (gamma_only) THEN
      CALL myzgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
               (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
                dvpsi, npwx )
   ELSE
      CALL myzgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
              dvpsi, npwx )
   END IF
ENDIF
!$acc end host_data
!
!$acc end data

IF (gamma_only) THEN
   !$acc exit data delete(ps_r)
   DEALLOCATE(ps_r)
ENDIF
DEALLOCATE(wwg)
DEALLOCATE(ps)
CALL stop_clock ('ortho')
RETURN
END SUBROUTINE orthogonalize_omega
