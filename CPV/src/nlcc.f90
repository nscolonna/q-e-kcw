!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   subroutine core_charge_ftr( tpre )
!=----------------------------------------------------------------------------=!
     !
     !  Compute the fourier transform of the core charge, from the radial
     !  mesh to the reciprocal space
     !
     use kinds,              ONLY : DP
     use ions_base,          ONLY : nsp
     use atom,               ONLY : rgrid
     use uspp,               ONLY : nlcc_any
     use uspp_param,         ONLY : upf
     use smallbox_gvec,      ONLY : ngb, gb
     use small_box,          ONLY : omegab, tpibab
     use pseudo_base,        ONLY : compute_rhocg
     use pseudopotential,    ONLY : tpstab, rhoc1_sp, rhocp_sp
     use cell_base,          ONLY : omega, tpiba2, tpiba
     USE splines,            ONLY : spline
     use gvect,              ONLY : gg, gstart
     USE core,               ONLY : rhocb, rhocg, drhocg
     USE fft_base,           ONLY: dfftp
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(IN) :: tpre
     !
     INTEGER :: is, ig
     REAL(DP) :: xg, cost1
     !
     !
     IF( .NOT. nlcc_any ) RETURN
     !
     IF( .NOT. ALLOCATED( rgrid ) ) &
        CALL errore( ' core_charge_ftr ', ' rgrid not allocated ', 1 )
     IF( .NOT. ALLOCATED( upf ) ) &
        CALL errore( ' core_charge_ftr ', ' upf not allocated ', 1 )
     !
     do is = 1, nsp
        !
        if( upf(is)%nlcc ) then
           !
              CALL compute_rhocg( rhocb(:,is), rhocb(:,is), rgrid(is)%r, &
                  rgrid(is)%rab, upf(is)%rho_atc(:), gb, omegab, tpibab**2, &
                  rgrid(is)%mesh, ngb, 0 )
              !
           IF( tpre ) THEN
              !
              IF( tpstab ) THEN
                 !
                 cost1 = 1.0d0/omega
                 !
                 IF( gstart == 2 ) THEN
                    rhocg (1,is) = rhoc1_sp(is)%y( 1 ) * cost1
                    drhocg(1,is) = 0.0d0
                 END IF
                 DO ig = gstart, SIZE( rhocg, 1 )
                    xg = SQRT( gg(ig) ) * tpiba
                    rhocg (ig,is) = spline( rhoc1_sp(is), xg ) * cost1
                    drhocg(ig,is) = spline( rhocp_sp(is), xg ) * cost1
                 END DO
                 !
              ELSE

                 CALL compute_rhocg( rhocg(:,is), drhocg(:,is), rgrid(is)%r, &
                                     rgrid(is)%rab, upf(is)%rho_atc(:), gg, &
                                     omega, tpiba2, rgrid(is)%mesh, dfftp%ngm, 1 )

              END IF
              !
           END IF
           !
        endif
        !
     end do

     return
   end subroutine core_charge_ftr



!-----------------------------------------------------------------------
      subroutine add_cc( rhoc, rhog, rhor )
!-----------------------------------------------------------------------
!
! add core correction to the charge density for exch-corr calculation
!
      USE kinds,              ONLY: DP
      use electrons_base,     only: nspin
      use control_flags,      only: iverbosity
      use io_global,          only: stdout
      use mp_global,          only: intra_bgrp_comm
      use cell_base,          only: omega
      USE mp,                 ONLY: mp_sum

      ! this isn't really needed, but if I remove it, ifc 7.1
      ! gives an "internal compiler error"
      use gvect, only: gstart
      USE fft_interfaces,     ONLY: fwfft
      USE fft_base,           ONLY: dfftp
      USE fft_helper_subroutines, ONLY: fftx_add_threed2oned_gamma
!
      implicit none
      !
      REAL(DP),    INTENT(IN)   :: rhoc( dfftp%nnr )
      REAL(DP),    INTENT(INOUT):: rhor( dfftp%nnr, nspin )
      COMPLEX(DP), INTENT(INOUT):: rhog( dfftp%ngm,  nspin )
      !
      COMPLEX(DP), ALLOCATABLE :: wrk1( : )
!
      integer :: ig, ir, iss, isup, isdw
      REAL(DP) :: rsum
      !
      IF( iverbosity > 1 ) THEN
         rsum = SUM( rhoc ) * omega / DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         CALL mp_sum( rsum, intra_bgrp_comm )
         WRITE( stdout, 10 ) rsum 
10       FORMAT( 3X, 'Core Charge = ', D14.6 )
      END IF
      !
      ! In r-space:
      !
      if ( nspin .eq. 1 ) then
         iss=1
         call daxpy(dfftp%nnr,1.d0,rhoc,1,rhor(1,iss),1)
      else
         isup=1
         isdw=2
         call daxpy(dfftp%nnr,0.5d0,rhoc,1,rhor(1,isup),1)
         call daxpy(dfftp%nnr,0.5d0,rhoc,1,rhor(1,isdw),1)
      end if 
      !
      ! rhoc(r) -> rhoc(g)  (wrk1 is used as work space)
      !
      allocate( wrk1( dfftp%nnr ) )

      wrk1(:) = rhoc(:)

      call fwfft('Rho',wrk1, dfftp )
      !
      ! In g-space:
      !
      if (nspin.eq.1) then
         CALL fftx_add_threed2oned_gamma( dfftp, wrk1, rhog(:,iss) )
      else
         wrk1 = wrk1 * 0.5d0
         CALL fftx_add_threed2oned_gamma( dfftp, wrk1, rhog(:,isup) )
         CALL fftx_add_threed2oned_gamma( dfftp, wrk1, rhog(:,isdw) )
      end if

      deallocate( wrk1 )
!
      return
      end subroutine add_cc


!
!-----------------------------------------------------------------------
      subroutine force_cc(irb,eigrb,vxc,fion1)
!-----------------------------------------------------------------------
!
!     core correction force: f = \int V_xc(r) (d rhoc(r)/d R_i) dr
!     same logic as in newd - uses box grid. For parallel execution:
!     the sum over node contributions is done in the calling routine
!
      USE kinds,             ONLY: DP, SP
      use electrons_base,    only: nspin
      use smallbox_gvec,     only: gxb, ngb
      use smallbox_subs,     only: fft_oned2box, boxdotgrid
      use cell_base,         only: omega
      use ions_base,         only: nsp, na, nat
      use small_box,         only: tpibab
      use uspp_param,        only: upf
      use core,              only: rhocb
      use fft_interfaces,    only: invfft
      use fft_base,          only: dfftb, dfftp
      use gvect,             only: gstart

      implicit none

! input
      integer, intent(in)        :: irb(3,nat)
      complex(dp), intent(in):: eigrb(ngb,nat)
      real(dp), intent(in)   :: vxc(dfftp%nnr,nspin)
! output
      real(dp), intent(inout):: fion1(3,nat)
! local
      integer :: iss, ix, ig, is, ia, nfft, isa
      real(dp) :: fac, res
      complex(dp) facg
      complex(sp), allocatable :: qv(:)
      complex(sp), allocatable :: fg1(:), fg2(:)
      real(dp), allocatable :: fcc(:,:)

#if defined(_OPENMP)
      INTEGER :: itid, mytid, ntids, omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
!
      call start_clock( 'forcecc' )

      fac = omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3*nspin)

!$omp parallel default(none) &      
!$omp          shared(nsp, na, ngb, eigrb, dfftb, irb, rhocb, &
!$omp                 gxb, nat, fac, upf, vxc, nspin, tpibab, fion1 ) &
!$omp          private(mytid, ntids, is, ia, nfft, ig, isa, qv, fg1, fg2, itid, res, ix, fcc, facg, iss )


      allocate( fcc( 3, nat ) )
      allocate( qv( dfftb%nnr ) )
      allocate( fg1( ngb ) )
      allocate( fg2( ngb ) )

      fcc(:,:) = 0.d0

      isa = 0

#if defined(_OPENMP)
      mytid = omp_get_thread_num()  ! take the thread ID
      ntids = omp_get_num_threads() ! take the number of threads
      itid  = 0
#endif

      do is = 1, nsp

         if( .not. upf(is)%nlcc ) then
            isa = isa + na(is) 
            cycle
         end if 

#if defined(__MPI)

         do ia = 1, na(is)
            nfft = 1
            if ( ( dfftb%np3( ia + isa ) <= 0 ) .OR. ( dfftb%np2( ia + isa) <= 0 ) ) cycle
#else
         !
         ! two fft's on two atoms at the same time (when possible)
         !
         do ia=1,na(is),2
            nfft=2
            if( ia .eq. na(is) ) nfft=1
#endif

#if defined(_OPENMP)
            IF ( mytid /= itid ) THEN
               itid = MOD( itid + 1, ntids )
               CYCLE
            ELSE
               itid = MOD( itid + 1, ntids )
            END IF
#endif

            do ix=1,3
               if (nfft.eq.2) then
                  do ig=1,ngb
                     facg = tpibab*CMPLX(0.d0,gxb(ix,ig),kind=DP)*rhocb(ig,is)
                     fg1(ig) = eigrb(ig,ia+isa  )*facg
                     fg2(ig) = eigrb(ig,ia+isa+1)*facg
                  end do
                  CALL fft_oned2box( qv, fg1, fg2 )
               else
                  do ig=1,ngb
                     facg = tpibab*CMPLX(0.d0,gxb(ix,ig),kind=DP)*rhocb(ig,is)
                     fg1(ig) = eigrb(ig,ia+isa)*facg
                  end do
                  CALL fft_oned2box( qv, fg1 )
               end if
!
               call invfft( qv, dfftb, ia+isa )
               !
               ! note that a factor 1/2 is hidden in fac if nspin=2
               !
               do iss=1,nspin
                  res = boxdotgrid(irb(:,ia  +isa),1,qv,vxc(:,iss))
                  fcc(ix,ia+isa) = fcc(ix,ia+isa) + fac * res
                  if (nfft.eq.2) then
                     res = boxdotgrid(irb(:,ia+1+isa),2,qv,vxc(:,iss))
                     fcc(ix,ia+1+isa) = fcc(ix,ia+1+isa) + fac * res 
                  end if
               end do
            end do
         end do

         isa = isa + na(is)

      end do

!
!$omp critical
      do ia = 1, nat
        fion1(:,ia) = fion1(:,ia) + fcc(:,ia)
      end do
!$omp end critical

      deallocate( qv )
      deallocate( fg1 )
      deallocate( fg2 )
      deallocate( fcc )

!$omp end parallel
!
      call stop_clock( 'forcecc' )

      return
      end subroutine force_cc


!
!-----------------------------------------------------------------------
      subroutine set_cc( irb, eigrb, rhoc )
!-----------------------------------------------------------------------
!
!     Calculate core charge contribution in real space, rhoc(r)
!     Same logic as for rhov: use box grid for core charges
! 
      use kinds, only: dp, sp
      use ions_base,         only: nsp, na, nat
      use uspp_param,        only: upf
      use smallbox_gvec,     only: ngb
      use smallbox_subs,     only: fft_oned2box, box2grid
      use control_flags,     only: iprint
      use core,              only: rhocb
      use fft_interfaces,    only: invfft
      use fft_base,          only: dfftb, dfftp

      implicit none
! input
      integer, intent(in)        :: irb(3,nat)
      complex(dp), intent(in):: eigrb(ngb,nat)
! output
      real(dp), intent(out)  :: rhoc(dfftp%nnr)
! local
      integer nfft, ig, is, ia, isa
      complex(sp), allocatable :: wrk1(:)
      complex(sp), allocatable :: qv(:)
      complex(sp), allocatable :: fg1(:), fg2(:)

#if defined(_OPENMP)
      INTEGER :: itid, mytid, ntids, omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif
!
      call start_clock( 'set_cc' )

      allocate( wrk1 ( dfftp%nnr ) )
      wrk1 (:) = (0.0, 0.0)
!
!$omp parallel default(none) &      
!$omp          shared(nsp, na, ngb, eigrb, dfftb, irb, rhocb, &
!$omp                 nat, upf, wrk1 ) &
!$omp          private(mytid, ntids, is, ia, nfft, ig, isa, qv, fg1, fg2, itid )

      allocate( qv ( dfftb%nnr ) )
      allocate( fg1 ( ngb ) )
      allocate( fg2 ( ngb ) )
!
      isa = 0

#if defined(_OPENMP)
      mytid = omp_get_thread_num()  ! take the thread ID
      ntids = omp_get_num_threads() ! take the number of threads
      itid  = 0
#endif

      do is = 1, nsp
         !
         if (.not.upf(is)%nlcc) then
            isa = isa + na(is)
            cycle
         end if
         !
#if defined(__MPI)
         do ia=1,na(is)
            nfft=1
            if ( ( dfftb%np3( ia + isa ) <= 0 ) .OR. ( dfftb%np2 ( ia + isa ) <= 0 ) ) cycle
#else
         !
         ! two ffts at the same time, on two atoms (if possible: nfft=2)
         !
         do ia=1,na(is),2
            nfft=2
            if( ia.eq.na(is) ) nfft=1
#endif

#if defined(_OPENMP)
            IF ( mytid /= itid ) THEN
               itid = MOD( itid + 1, ntids )
               CYCLE
            ELSE
               itid = MOD( itid + 1, ntids )
            END IF
#endif

            if(nfft.eq.2)then
               fg1 = eigrb(1:ngb,ia  +isa)*rhocb(1:ngb,is)
               fg2 = eigrb(1:ngb,ia+1+isa)*rhocb(1:ngb,is)
               CALL fft_oned2box( qv, fg1, fg2 )
            else
               fg1 = eigrb(1:ngb,ia  +isa)*rhocb(1:ngb,is)
               CALL fft_oned2box( qv, fg1 )
            endif
!
            call invfft( qv, dfftb, isa+ia )
!
            call box2grid(irb(:,ia+isa),1,qv,wrk1)
            if (nfft.eq.2) call box2grid(irb(:,ia+1+isa),2,qv,wrk1)
!
         end do
         isa = isa + na(is)
      end do
!
      deallocate( qv  )
      deallocate( fg1  )
      deallocate( fg2  )

!$omp end parallel

      call copy( rhoc, wrk1 )

      deallocate( wrk1 )
!
      call stop_clock( 'set_cc' )
!
      return

   contains

      SUBROUTINE copy(rhoc,wsp,wdp)
         real(dp) :: rhoc(:)
         complex(dp), optional :: wdp(:)
         complex(sp), optional :: wsp(:)
         if(present(wdp)) then 
            call dcopy( dfftp%nnr, wdp, 2, rhoc, 1 )
         else if(present(wsp)) then
            call scopy( dfftp%nnr, wsp, 2, rhoc, 1 )
         endif
         RETURN
      END SUBROUTINE 

   end subroutine set_cc

