SUBROUTINE kcw_kpoint_grid()
!
! this routine calculates xq_ibz, nqstot_ibz, wq_ibz, ibz2fbz and fbz2ibz for 
! all the wannier functions. 
!
  USE kinds,              ONLY : DP  
  USE io_global,          ONLY : stdout
  USE control_kcw,        ONLY : num_wann,nsym_w_q, s_w
  USE control_kcw,        ONLY : nqstot_ibz, xq_ibz, wq_ibz
  USE control_kcw,        ONLY : mp1, mp2, mp3
  USE control_kcw,        ONLY : kcw_iverbosity
  USE klist,              ONLY : xk
  USE control_kcw,        ONLY : ibz2fbz, fbz2ibz
  USE cell_base,          ONLY : bg
  USE symm_base,          ONLY : time_reversal, t_rev, sr
  USE control_kcw,        ONLY:   sym_w2sym, nsym_w_q
  USE control_kcw,        ONLY : nqstar, qstar_iq, qstar_isym, qstar_Gvec, nqstot
  !
  implicit none
  !
  INTEGER                :: iwann, iq, isym, iRxq
  ! indices
  INTEGER                :: nsym_w_iwann
  ! number of symmetries respected by current wannier function
  INTEGER, ALLOCATABLE   :: s_w_iwann(:,:,:)
  ! subset of s, respected by current wannier function iwann
  REAL(DP)               :: xq(3)
  !current q point
  REAL(DP)               :: Gvector(3), Gvector_cryst(3)
  ! G vector connecting xq with xk
  REAL(DP), ALLOCATABLE  :: xq_ibz_iwann(:,:), wq_ibz_iwann(:)
  !
  ALLOCATE( nqstot_ibz( num_wann ) )
  ALLOCATE( xq_ibz( 3, mp1*mp2*mp3, num_wann ) )
  ALLOCATE( wq_ibz( mp1*mp2*mp3, num_wann ) )
  ALLOCATE( xq_ibz_iwann( 3, mp1*mp2*mp3 ) )
  ALLOCATE( wq_ibz_iwann( mp1*mp2*mp3 ) )
  ALLOCATE( s_w_iwann (3,3,48) )
  ALLOCATE( ibz2fbz( mp1*mp2*mp3, num_wann))
  ALLOCATE( fbz2ibz( mp1*mp2*mp3, num_wann))  
  !
  !WRITE(stdout,*) "Finding symmetries of wannier functions......"
  WRITE(stdout,'(/, 5X, "SYM : Finding the IBZ")') 
  DO iwann = 1, num_wann
    nsym_w_iwann = nsym_w_q(iwann)
    s_w_iwann(:,:,:) = s_w(:,:,:,iwann)
    CALL kpoint_grid ( nsym_w_iwann, time_reversal, .false., s_w_iwann, 0, bg, &
                       mp1*mp2*mp3, 0,0,0, mp1,mp2,mp3, &
                       nqstot_ibz(iwann), xq_ibz_iwann, wq_ibz_iwann )
    WRITE(stdout,'(7X, "iwann =", I5, 3X, "nqstot_ibz =", I5, 3X)') iwann, nqstot_ibz(iwann) 
    !
    !in general, kpoint_grid gives points in other BZ, we find the index
    !in our original FBZ (xk) so that we are able to map xq in xk
    !
    fbz2ibz(:, iwann) = -1
    DO iq = 1, nqstot_ibz(iwann)
      xq(:) = xq_ibz_iwann(:, iq)
      CALL find_kpoint(xq, ibz2fbz(iq,iwann), Gvector, Gvector_cryst)
      xq_ibz(:, iq, iwann) = xk( :, ibz2fbz(iq, iwann) )
      wq_ibz( iq, iwann ) = wq_ibz_iwann (iq)
      fbz2ibz( ibz2fbz(iq, iwann), iwann ) = iq
    END DO!iq
    !
    IF (kcw_iverbosity .gt. 1) THEN
      WRITE(stdout,'(9X, "xq(iq="i3, x, ") = ", 3F8.4, 3X, "wq = ", F10.6, 3X, "iq_FBZ = ", I3)') & 
              (iq, xq_ibz(:, iq, iwann), wq_ibz( iq, iwann ), ibz2fbz(iq,iwann), iq=1,nqstot_ibz(iwann))
      WRITE(stdout, *)
    ENDIF 
  END DO!iwann
  !
  DEALLOCATE( xq_ibz_iwann )
  DEALLOCATE( wq_ibz_iwann )
  DEALLOCATE( s_w_iwann )
  !
END SUBROUTINE



SUBROUTINE kcw_starq()
  USE kinds,              ONLY : DP  
  USE control_kcw,        ONLY : nqstar, qstar_iq, qstar_isym, qstar_Gvec, nqstot, xq_ibz
  USE control_kcw,        ONLY : nqstot_ibz,ibz2fbz, nsym_w_q, sym_w2sym, num_wann
  USE symm_base,          ONLY : time_reversal, t_rev, sr, nsym
  !
  ! generate the star of q for each q point
  !
  INTEGER :: iq_ibz, iq, iwann
  REAL(DP) :: xq(3)
  REAL(DP)               :: Gvector(3), Gvector_cryst(3)
  !
  WRITE(*,*) "starq", "nqstot_ibz:", nqstot_ibz, "xq_ibz", xq_ibz
  ALLOCATE ( nqstar(nqstot, num_wann) )
  ALLOCATE ( qstar_iq(nqstot, nqstot, num_wann) )
  ALLOCATE ( qstar_isym(nqstot, nqstot, num_wann) )
  ALLOCATE ( qstar_Gvec(3, nqstot, nqstot, num_wann) )
  !
  nqstar(:,:) = 0
  DO iwann = 1, num_wann
    CALL reset_symmetry_op(iwann)
    WRITE(*,*) "iwann: ", iwann, "nqstot_ibz(iwann)", nqstot_ibz(iwann)
    DO iq_ibz = 1, nqstot_ibz(iwann)
      iq = ibz2fbz(iq_ibz,iwann)
      xq(:) = xq_ibz(:, iq_ibz, iwann)
      WRITE(*,*) "xq_ibz", xq
      WRITE(*,*) "nsym", nsym
      DO isym = 1, nsym
        WRITE(*,*) "isym", isym
        CALL rotate_xk(xq, iRxq, Gvector, Gvector_cryst, &
                       sr(:,:,isym), t_rev(isym))
                       WRITE(*,*) "iRxq", iRxq
        ! check if this point is already in the list of the star
        IF( .NOT. ( ANY(qstar_iq(:, iq, iwann) .EQ. iRxq ) ) ) THEN
          WRITE(*,*) "nqstar"
          nqstar(iq, iwann) = nqstar(iq, iwann) + 1
          WRITE(*,*) "qstar_iq"
          qstar_iq(nqstar(iq, iwann), iq, iwann) = iRxq
          WRITE(*,*) "qstar_isym"
          qstar_isym(nqstar(iq, iwann), iq, iwann) = isym
          WRITE(*,*) "qstar_Gvec"
          qstar_Gvec(:, nqstar(iq, iwann), iq, iwann) = Gvector(:)
        END IF
      END DO!isym
    END DO!iq
    WRITE(*,*) "nqstar:", nqstar
    WRITE(*,*) "qstar_iq", qstar_iq
  END DO!iwann
END SUBROUTINE
