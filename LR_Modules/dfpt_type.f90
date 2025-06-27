!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
MODULE dfpt_type
   !
   USE kinds, ONLY : DP
   !
   IMPLICIT NONE
   !
   ! Data that describes linear response quantities.
   !
   TYPE :: dfpt_data_type
      INTEGER :: npert
      !! Number of perturbations
      COMPLEX(DP), ALLOCATABLE :: drhos(:, :, :)
      !! Change of the charge density (smooth part only, dffts)
      !! If norm conserving PP and doublegrid = .FALSE., identical to drhop.
      !! If USPP or PAW, does not include the augmentation charges. (drhop does.)
      COMPLEX(DP), ALLOCATABLE :: drhop(:, :, :)
      !! Change of the charge density (smooth and hard parts, dfftp)
      COMPLEX(DP), ALLOCATABLE :: dvscfs(:, :, :)
      !! Change of the scf potential (smooth part only, dffts)
      !! If doublegrid = .FALSE., identical to dvscfp, even for USPP or PAW.
      COMPLEX(DP), ALLOCATABLE :: dvscfp(:, :, :)
      !! Change of the scf potential (smooth and hard parts, dfftp)
      COMPLEX(DP), ALLOCATABLE :: dbecsum(:, :, :, :)
      !! Change of becsum. Only used for USPP and PAW.
      !
      ! Variables allocated and used only for phonon perturbations
      !
      COMPLEX(DP), ALLOCATABLE :: drhoc(:, :)
      !! Change of core charge due to nonlinear core correction.
      !! Size (dfftp%nnr, npert)
      COMPLEX(DP), ALLOCATABLE :: drhop_pulay(:, :, :)
      !! Pulay correction to drhop due to augmentation charge. Used for USPP or PAW.
      !! Size (dfftp%nnr, nspin_mag, npert)
      COMPLEX(DP), ALLOCATABLE :: dbecsum_pulay(:, :, :, :)
      !! Pulay correction to dbecsum due to augmentation charge. Used for PAW only.
      !! Size ((nhm * (nhm + 1))/2, nat, nspin_mag, npert)
   END TYPE dfpt_data_type
   !
   CONTAINS
   !
   !---------------------------------------------------------------------------
   SUBROUTINE dfpt_dvscfp_to_dvscfs(dfpt_data)
      USE fft_base,              ONLY : dffts, dfftp
      USE fft_interfaces,        ONLY : fft_interpolate
      USE gvecs,                 ONLY : doublegrid
      USE noncollin_module,      ONLY : nspin_mag
      !
      IMPLICIT NONE
      !
      TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
      !
      INTEGER :: ipert, is
      !
      IF (doublegrid) THEN
         DO ipert = 1, dfpt_data%npert
            DO is = 1, nspin_mag
               CALL fft_interpolate(dfftp, dfpt_data%dvscfp(:, is, ipert), &
                                    dffts, dfpt_data%dvscfs(:, is, ipert))
            ENDDO
         ENDDO
      ELSE
         ! If doublegrid is false, dvscfs and dvscfp are the same.
         CALL zcopy(dfftp%nnr*nspin_mag*dfpt_data%npert, dfpt_data%dvscfp, 1, dfpt_data%dvscfs, 1)
      ENDIF
      !
   END SUBROUTINE dfpt_dvscfp_to_dvscfs
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   SUBROUTINE allocate_dfpt_data(dfpt_data, npert)
      !
      USE fft_base,             ONLY : dfftp, dffts
      USE noncollin_module,     ONLY : nspin_mag
      USE uspp_param,           ONLY : nhm
      USE ions_base,            ONLY : nat
      USE uspp,                 ONLY : okvan
      !
      IMPLICIT NONE
      !
      TYPE(dfpt_data_type), INTENT(OUT) :: dfpt_data
      INTEGER, INTENT(IN) :: npert
      !
      dfpt_data%npert = npert
      !
      ALLOCATE(dfpt_data%drhos(dffts%nnr, nspin_mag, npert))
      ALLOCATE(dfpt_data%drhop(dfftp%nnr, nspin_mag, npert))
      ALLOCATE(dfpt_data%dvscfs(dffts%nnr, nspin_mag, npert))
      ALLOCATE(dfpt_data%dvscfp(dfftp%nnr, nspin_mag, npert))
      dfpt_data%drhos = (0.d0, 0.d0)
      dfpt_data%drhop = (0.d0, 0.d0)
      dfpt_data%dvscfs = (0.d0, 0.d0)
      dfpt_data%dvscfp = (0.d0, 0.d0)
      !
      IF (okvan) THEN
         ALLOCATE(dfpt_data%dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag, npert))
         dfpt_data%dbecsum = (0.d0, 0.d0)
      ELSE
         ! If okvan is false, dbecsum is not used.
         ! TODO: Do not allocate it (need to change the code that uses it).
         ALLOCATE(dfpt_data%dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag, npert))
      ENDIF
      !
   END SUBROUTINE allocate_dfpt_data
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   SUBROUTINE deallocate_dfpt_data(dfpt_data)
      !
      TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
      !
      IF (ALLOCATED(dfpt_data%drhos)) DEALLOCATE(dfpt_data%drhos)
      IF (ALLOCATED(dfpt_data%drhop)) DEALLOCATE(dfpt_data%drhop)
      IF (ALLOCATED(dfpt_data%dvscfp)) DEALLOCATE(dfpt_data%dvscfp)
      IF (ALLOCATED(dfpt_data%dvscfs)) DEALLOCATE(dfpt_data%dvscfs)
      IF (ALLOCATED(dfpt_data%dbecsum)) DEALLOCATE(dfpt_data%dbecsum)
      IF (ALLOCATED(dfpt_data%drhoc)) DEALLOCATE(dfpt_data%drhoc)
      IF (ALLOCATED(dfpt_data%drhop_pulay)) DEALLOCATE(dfpt_data%drhop_pulay)
      IF (ALLOCATED(dfpt_data%dbecsum_pulay)) DEALLOCATE(dfpt_data%dbecsum_pulay)
      !
   END SUBROUTINE deallocate_dfpt_data
   !---------------------------------------------------------------------------
END MODULE dfpt_type
