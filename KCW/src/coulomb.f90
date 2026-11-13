MODULE coulomb
  !
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER       :: nq1=1, nq2=1, nq3=1
  REAL(DP)      :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP)      :: eps  = 1.d-6
  REAL(DP)      :: exxdiv     = 0._dp
  REAL(DP)      :: exxdiv_eps = 0._dp
  CHARACTER(32) :: exxdiv_treatment  = ' '
  LOGICAL       :: use_regularization = .FALSE.
  ! ... x_gamma_extrapolation
  LOGICAL       :: x_gamma_extrapolation =.FALSE.
  LOGICAl       :: on_double_grid =.FALSE.
  REAL(DP)      :: grid_factor = 1.d0 !8.d0/7.d0
  !
  REAL(DP)      :: yukawa = 0._dp
  REAL(DP)      :: erfc_scrlen = 0._dp
  REAL(DP)      :: erf_scrlen = 0._dp
  REAL(DP)      :: gau_scrlen = 0.d0
  REAL(DP)      :: eps_mat(3,3)
  !
  CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE setup_coulomb ( )
    !--------------------------------------------------------------------
    !
    !	
    !!  This routine initialize the coulomb kernel according to the 
    !!  Gygi-Balderschi scheme. This is compatible with peridic 
    !!  systems and finite sytems calculation. For finite system 
    !!  one can also (should) use assume_isolated (see kcw_readin) 
    !!  Adapted from exx_divergence () inside PW. Generalized to the
    !!  screened Coulomb eps^-1 * vc adpting the madelung constant 
    !!  formula in Rurali and Cartoixa Nano Letters 9. 975 (2009)
    !
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kcw,          ONLY : mp1, mp2, mp3, l_vcut, eps_inf, calculation
    USE martyna_tuckerman,    ONLY : do_comp_mt
    USE cell_base,            ONLY : omega
    USE Coul_cut_2D,          ONLY : do_cutoff_2D
    !
    IMPLICIT NONE 
    !
    LOGICAL :: exst, skip_eps
    INTEGER :: nqs
    !
    CALL start_clock( 'Coulomb setup' )
    !
    nq1=mp1; nq2=mp2; nq3=mp3
    nqs=nq1*nq2*nq3
    !
    eps_mat = 0.D0
    eps_mat (1,1) = 1.D0; eps_mat (2,2) = 1.D0; eps_mat (3,3) = 1.D0
    !
    yukawa = 0._dp
    erfc_scrlen = 0._dp
    erf_scrlen = 0._dp
    gau_scrlen = 0.d0
    exxdiv_treatment = 'none'
    x_gamma_extrapolation = .false.
    ! NOTABENE: The Screened contribution does not work with x_gamma_extrapolation
    !           Leave x_gamma_extrapolation to FALSE unless you know what you are doing
    !
    skip_eps = .FALSE. 
    IF (l_vcut) THEN
      !
      use_regularization = .true. 
      exxdiv_treatment = 'gb' ! For now only GB scheme
      !
      ! Set the dielectric tensor. Priority is given in this order
      !   1) if a file eps.dat exist, read and set up the dielectric tensor
      INQUIRE( file="eps.dat", exist=exst )
      IF ( exst) THEN  
         !
         WRITE (stdout,'(/,5X, "INFO: Dielectric tensor read from file eps.dat")')  
         CALL read_eps (eps_mat)
         eps_inf = (eps_mat(1,1)+eps_mat(2,2)+eps_mat(3,3))/3.D0
         WRITE (stdout,'(  5X, "INFO: average macroscopic eps", 1F12.6)')  eps_inf
         !
         ! 2) the average epsilon_inf is given in input
      ELSE IF ( ABS(eps_inf -1.D0) .gt. 1E-6) THEN
         WRITE (stdout,'(/,5X, "INFO: average macroscopic eps from input", 1F12.6)')  eps_inf
         eps_mat (1,1) = eps_inf; eps_mat (2,2) = eps_inf; eps_mat (3,3) = eps_inf
         !
         ! 3) No information about the dielectric properties: set eps=I
      ELSE
         ! With the 2D cutoff the screened q+G=0 term cannot simply be dropped: it
         ! contains the (unscreened) charged-sheet term that cancels a vacuum-dependent
         ! artefact, and its in-plane part depends strongly on alpha_2D. Require eps.dat
         ! whenever the screened term is used (screen) or printed (cc).
         IF (do_cutoff_2D .AND. (calculation == 'screen' .OR. calculation == 'cc')) &
            CALL errore('setup_coulomb', '2D cutoff with l_vcut=.true.: provide the supercell &
                        &dielectric tensor in eps.dat (its in-plane block gives alpha_2D)', 1)
         WRITE (stdout,'(/,5X, "INFO: Dielectric tensor = NOT SPECIFIED")')
         IF (.NOT. do_cutoff_2D) &
            WRITE (stdout,'(  5X, "      NO Correction for the Screened Coulomb")')
         skip_eps = .TRUE.
         !
      ENDIF
      !
    ENDIF
    !
    CALL divergence ()
    IF (skip_eps) exxdiv_eps = 0.D0
    !
    ! With the 2D cutoff both the bare and the screened q+G=0 terms are recomputed
    ! (charged-sheet term + 2D Gygi-Baldereschi, see divergence_2d). The 3D GB above is
    ! the correction for the 3D-periodic kernel (a stack of layers) and is inconsistent
    ! with the cut-off kernel; a constant eps tensor is not valid in 2D.
    IF (do_cutoff_2D .AND. use_regularization) CALL divergence_2d ( skip_eps )
    ! Compute the divergence and set the q+G=0 term (exxdiv and exxdiv_eps, 
    ! to be used to correct the bare and screened coulomb. USed in bare_pot
    ! and screen_coeff. 
    !
    IF (.NOT. do_comp_mt) THEN
      !
      WRITE (stdout,'(/,5X, "INFO: Coulomb q+G=0 treatment:")') 
      WRITE (stdout,'(  5X, "INFO: Divergence         ", 3x, 1A8)'     )  exxdiv_treatment
      WRITE( stdout,'(  5X, "INFO: q-grid dimension   ", 3x, 3I4)'     )  nq1, nq2, nq3
      WRITE( stdout,'(  5X, "INFO: cell volume        ", 3x, 1F20.12)' )  omega
      WRITE (stdout,'(  5X, "INFO: Gamma Extrapolation", 3x, 1L5 )'    )  x_gamma_extrapolation
      !
      IF ( x_gamma_extrapolation ) THEN
         !
         WRITE( stdout, '(5X, "INFO: extrapolation q->0 dealt with 8/7 -1/7 trick ")')
         grid_factor = 8.d0 / 7.d0
         !
      ELSE
         !
         WRITE( stdout, '(5X, "INFO: extrapolation q->0 term not estimated")' )
         grid_factor = 1.d0
         !
      ENDIF
      !
      WRITE (stdout,'(  5X,    "INFO: Bare Coulomb q+G=0     ", 3x, 1ES15.5 )')  exxdiv
      IF ( (calculation == "screen" .AND. l_vcut) .OR. calculation == 'cc' ) THEN
         !
         IF (.NOT. skip_eps) THEN 
           WRITE (stdout,'(  5X, "INFO: Epsilon infinity       ", 3x, 1F12.6  )')  eps_inf
           WRITE (stdout,'(  5X, "INFO: Dielectric tensor      ", 3x, 3F12.6  )') (eps_mat(:,1))
           WRITE (stdout,'(  5X, "                             ", 3x, 3F12.6  )') (eps_mat(:,2))
           WRITE (stdout,'(  5X, "                             ", 3x, 3F12.6  )') (eps_mat(:,3))
         ENDIF
         WRITE (stdout,'(  5X, "INFO: Screened Coulomb q+G=0 ", 3x, 1F20.12 )')  exxdiv_eps
         IF (.NOT. skip_eps) THEN 
            WRITE (stdout,'(  5X, "      Isotropic estimate     ", 3x, 1F20.12 )')  exxdiv/eps_inf
         ELSE
            WRITE (stdout,'(  5X, "      Isotropic estimate     ", 3x, 1F20.12 )')  0.D0
         ENDIF
         WRITE (stdout,'(  5X, "INFO: uPi(q+G=0) estimation  ", 3X, 1F20.12 )') -exxdiv/omega/nqs
         WRITE (stdout,'(  5X, "INFO: rPi(q+G=0) estimation  ", 3X, 1F20.12 )') -(exxdiv_eps)/omega/nqs
         !
      ENDIF
      !
    ELSE 
      !
      WRITE (stdout,'(/,5X, "INFO: Coulomb kernel treated according to MT scheme")') 
      !
    ENDIF 
    !
    CALL stop_clock( 'Coulomb setup' )
    !  
    RETURN
    !
  END SUBROUTINE setup_coulomb
  !
  !----------------------------------------------------------------------
  SUBROUTINE g2_convolution( ngm, g, xk, xkq, fac )
    !--------------------------------------------------------------------
    !! This routine calculates the 1/|r-r'| part of the exact exchange
    !! expression in reciprocal space (the G^-2 factor).
    !! It then regularizes it according to the specified recipe.
    !
    USE kinds,       ONLY : DP
    USE cell_base,   ONLY : at, tpiba2, tpiba, alat
    USE constants,   ONLY : fpi, e2, pi
    USE Coul_cut_2D, ONLY : do_cutoff_2D
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: ngm
    !! Number of G vectors
    REAL(DP), INTENT(IN) :: g(3,ngm)
    !! Cartesian components of G vectors
    REAL(DP), INTENT(IN) :: xk(3)
    !! current k vector
    REAL(DP), INTENT(IN) :: xkq(3)
    !! current q vector
    REAL(DP), INTENT(INOUT) :: fac(ngm)
    !! Calculated convolution
    !
    ! ... local variables
    !
    INTEGER :: ig !Counters
    REAL(DP) :: q(3), qq, x
    REAL(DP) :: grid_factor_track(ngm), qq_track(ngm)
    REAL(DP) :: nqhalf_dble(3)
    LOGICAL :: odg(3)
    REAL(DP) :: cut2d(ngm), lz
    ! ... 2D cutoff factor 1 - exp(-|q+G|_par lz) cos((q+G)_z lz), lz = c/2
    !     (PRB 96, 075448, Eq. 24). Set to 1 if the 2D cutoff is not active.
    !
    IF (do_cutoff_2D) THEN
       lz = 0.5_DP*at(3,3)*alat
!$omp parallel do default(shared), private(ig,q)
       DO ig = 1, ngm
          q(:) = xk(:) - xkq(:) + g(:,ig)
          cut2d(ig) = 1._DP - EXP(-SQRT(q(1)**2+q(2)**2)*tpiba*lz) * COS(q(3)*tpiba*lz)
       ENDDO
!$omp end parallel do
    ELSE
       cut2d(:) = 1._DP
    ENDIF
    !
    ! Now the Coulomb potential that are computed on the fly
    !
    nqhalf_dble(1:3) = (/ DBLE(nq1)*0.5_DP, DBLE(nq2)*0.5_DP, DBLE(nq3)*0.5_DP /)
    !
    ! Set the grid_factor_track and qq_track
    !
    IF ( x_gamma_extrapolation ) THEN
!$omp parallel do default(shared), private(ig,q,x,odg)
       DO ig = 1, ngm
          q(:)= xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = SUM(q(:)**2) * tpiba2
          x = (q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nqhalf_dble(1)
          odg(1) = ABS(x-NINT(x)) < eps
          x = (q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nqhalf_dble(2)
          odg(2) = ABS(x-NINT(x)) < eps
          x = (q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nqhalf_dble(3)
          odg(3) = ABS(x-NINT(x)) < eps
          IF( ALL( odg(:) ) ) THEN
             grid_factor_track(ig) = 0._DP ! on double grid
          ELSE
             grid_factor_track(ig) = grid_factor ! not on double grid
          ENDIF
       ENDDO
!$omp end parallel do
    ELSE
!$omp parallel do default(shared), private(ig,q)
       DO ig = 1, ngm
          q(:) = xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = SUM(q(:)**2) * tpiba2
       ENDDO
       grid_factor_track = 1._DP
       !
    ENDIF
    ! ... The big loop
    !
!$omp parallel do default(shared), private(ig,qq)
    DO ig = 1, ngm
      !
      qq = qq_track(ig)
      !
      IF(gau_scrlen > 0) THEN
         fac(ig) = e2*((pi/gau_scrlen)**(1.5_DP))*EXP(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
         !
      ELSEIF (qq > eps_qdiv) THEN
         !
         IF ( erfc_scrlen > 0  ) THEN
            fac(ig) = e2*fpi/qq*(1._DP-EXP(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
         ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig) = e2*fpi/qq*(EXP(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
         ELSE
            fac(ig) = e2*fpi/( qq + yukawa ) * grid_factor_track(ig) * cut2d(ig) ! as HARTREE (2D cutoff if active)
         ENDIF
         !
      ELSE
         !
         fac(ig) = - exxdiv ! or rather something ELSE (see F.Gygi)
         !
         IF (yukawa>0._DP .AND. .NOT.x_gamma_extrapolation) fac(ig) = fac(ig) + &
                                                            e2*fpi/( qq + yukawa )
         !
         IF (erfc_scrlen>0._DP .AND. .NOT.x_gamma_extrapolation) fac(ig) = fac(ig) + &
                                                                  e2*pi/(erfc_scrlen**2)
         !
      ENDIF
      !
    ENDDO
!$omp end parallel do
    !
  END SUBROUTINE g2_convolution
  !
  !-----------------------------------------------------------------------
  SUBROUTINE divergence()
      !-----------------------------------------------------------------------
    !
    ! ...  This function calculates the G=0 term of the Coulomb potential
    ! ...  in presence  of a unitary charge (q=1). Used for counter-charge
    ! ...  corrections.
    !
    USE kinds,              ONLY : DP
    USE constants,          ONLY : fpi, e2
    USE cell_base,          ONLY : bg, at, alat, omega
    USE gvect,              ONLY : ngm, g
    USE gvecw,              ONLY : gcutw
    USE control_flags,      ONLY : gamma_only
    USE mp_global,          ONLY : intra_pool_comm
    USE mp,                 ONLY : mp_sum
    USE control_kcw,        ONLY : eps_inf
    !
    !
    IMPLICIT NONE
    !
    INTEGER  :: nqs                   ! number of points in the q-grid
    LOGICAL  :: on_double_grid=.false.
    INTEGER  :: iq1,iq2,iq3, ig
    REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)
    INTEGER  :: nqq, iq
    REAL(DP) :: aa, dq
    !
    REAL(DP) :: q_eps_q
    REAL(DP) :: aa_eps, div_eps, det_eps, alpha_eps, dq_eps
    !
    CALL start_clock( 'exx_div' )
    !
    IF ( .NOT. use_regularization ) THEN
       exxdiv     = 0._DP
       exxdiv_eps = 0._DP
       RETURN
    ENDIF
    !
    nqs = nq1 * nq2 * nq3
    tpiba2 = ( fpi / 2.d0 / alat ) ** 2
    alpha  = 10.d0 / gcutw
    alpha_eps = alpha/eps_inf
    !
    dq1= 1.d0/DBLE(nq1)
    dq2= 1.d0/DBLE(nq2) 
    dq3= 1.d0/DBLE(nq3)
    !
    div = 0.d0
    div_eps = 0.d0
    !
    DO iq1 = 1, nq1
      DO iq2 = 1, nq2
        DO iq3 = 1, nq3
          !
          xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                  bg(:,2) * (iq2-1) * dq2 + &
                  bg(:,3) * (iq3-1) * dq3
          !
          DO ig = 1, ngm
            !
            q(1) = xq(1) + g(1,ig)
            q(2) = xq(2) + g(2,ig)
            q(3) = xq(3) + g(3,ig)
            qq = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
            q_eps_q = DOT_PRODUCT( q(:), MATMUL( eps_mat(:,:), q(:) ) )
            !
            IF ( x_gamma_extrapolation ) THEN
              !
              on_double_grid = .true.
              x = 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
              on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
              x = 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
              on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
              x = 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
              on_double_grid = on_double_grid .AND. (abs(x-nint(x))<eps)
              !
            ENDIF
            !
            IF ( .NOT. on_double_grid ) THEN
              IF ( qq > 1.d-8 ) THEN
                div = div + EXP( - alpha * qq ) / ( qq + yukawa / tpiba2 ) &
                                                    * grid_factor
                div_eps = div_eps + EXP( - alpha_eps * q_eps_q ) / ( q_eps_q + yukawa / tpiba2 ) &
                                                    * grid_factor
              ENDIF
            ENDIF
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum( div, intra_pool_comm )
    CALL mp_sum( div_eps, intra_pool_comm )
    !
    IF ( gamma_only ) div = 2.d0 * div
    IF ( gamma_only ) div_eps = 2.d0 * div_eps
    !
    IF ( .NOT. x_gamma_extrapolation ) THEN
      IF ( yukawa < 1.d-8) THEN
        div = div - alpha
        div_eps = div_eps - alpha_eps
      ELSE
        div = div + tpiba2 / yukawa
        div_eps = div_eps + tpiba2 / yukawa
      ENDIF
    ENDIF
    !
    div = div * e2 * fpi / tpiba2 / nqs
    div_eps = div_eps * e2 * fpi / tpiba2 / nqs
    !
    alpha = alpha / tpiba2
    alpha_eps = alpha_eps / tpiba2
    !
    nqq = 100000
    dq = 5.0d0 / SQRT( alpha ) / nqq
    dq_eps = 5.0d0 / SQRT( alpha_eps ) / nqq
    aa = 0.d0
    aa_eps = 0.d0
    !
    DO iq = 0, nqq
      !
      q_ = dq * ( iq + 0.5d0 )
      qq = q_ * q_
      aa = aa - EXP( - alpha * qq ) * yukawa / ( qq + yukawa ) * dq
      q_ = dq_eps * ( iq + 0.5d0 )
      qq = q_ * q_
      aa_eps = aa_eps - EXP( - alpha_eps * q_eps_q ) * yukawa / ( q_eps_q + yukawa ) * dq_eps
      !
    ENDDO
    !
    det_eps = eps_mat(1,1)*(eps_mat(2,2)*eps_mat(3,3)-eps_mat(2,3)*eps_mat(3,2)) & 
             -eps_mat(2,2)*(eps_mat(2,1)*eps_mat(3,3)-eps_mat(2,3)*eps_mat(3,1)) &
             +eps_mat(3,3)*(eps_mat(2,1)*eps_mat(3,2)-eps_mat(2,2)*eps_mat(3,1)) 
    !
    aa = aa * 8.d0 / fpi
    aa_eps = aa_eps * 8.d0 / fpi
    aa = aa + 1.d0 / SQRT( alpha * 0.25d0 * fpi )
    aa_eps = aa_eps + 1.d0 / SQRT( alpha_eps * 0.25d0 * fpi * det_eps)
    !  
    div = div - e2 * omega * aa
    div_eps = div_eps - e2 * omega * aa_eps
    !
    exxdiv = div * nqs
    exxdiv_eps = div_eps * nqs
    !
    CALL stop_clock( 'exx_div' )
    !
    RETURN
    !
  END SUBROUTINE divergence
  !
  !----------------------------------------------------------------------
  SUBROUTINE divergence_2d( skip_eps )
    !--------------------------------------------------------------------
    !
    !! Bare and screened q+G=0 terms for a 2D system with the 2D Coulomb cutoff
    !! (lz = c/2, mesh nq1 x nq2 x 1, in-plane cell area A, Omega = 2 lz A).
    !!
    !! 1) Charged-sheet term. At q=0 the code drops G=0 but keeps G=(0,0,Gz). With the
    !!    cutoff kernel, these Gz/=0 terms of a sheet of total charge 1 per supercell
    !!    contain the shape-independent piece  pi e2 lz / A  (the background spread over
    !!    the whole cell, vacuum included). It is cancelled by -2pi e2 lz^2 in fac(0).
    !!    It is NOT screened: at q=0 the net sheet charge can only move along z.
    !! 2) In-plane 1/q singularity: 2D Gygi-Baldereschi. The q->0 limit of the full Gz
    !!    sum is the isolated sheet 2pi e2/(q A), i.e. a G=0 kernel 4pi e2 lz/q. Choose
    !!    the q=0 value F0 so that the FULL grid sum of an auxiliary function with this
    !!    singularity reproduces its integral:
    !!        F(q)   = 4pi e2 lz exp(-beta q^2) / q                       (bare)
    !!        F_s(q) = 4pi e2 lz exp(-beta q^2) / (q (1 + r0(qhat) q))    (screened)
    !!        F0     = Nq A/(2pi)^2 int d^2q F  -  sum'_{q in grid, G_par} F(q+G_par)
    !!    The remainder (true kernel - auxiliary) is smooth, so its grid sum converges
    !!    fast. Averaging 1/q over the central cell only (previous version) leaves an
    !!    O(1/nq) error.
    !!    eps(q) = 1 + r0 q (strict 2D, Keldysh), r0 = 2pi qhat.alpha_2D.qhat, with
    !!    alpha_2D(i,j) = (eps_mat(i,j) - delta_ij) c/(4pi), i,j=x,y, from the in-plane
    !!    block of the SUPERCELL dielectric tensor (eps.dat).
    !! KCW convention: fac(q+G=0) = -exxdiv, fac_eps(q+G=0) = -exxdiv_eps, so
    !!        -exxdiv = F0 - 2pi e2 lz^2,     -exxdiv_eps = F0_s - 2pi e2 lz^2
    !! beta = 10/ecutwfc (bohr^2), the same Gaussian width as the 3D GB in divergence().
    !
    USE kinds,              ONLY : DP
    USE constants,          ONLY : fpi, pi, e2, tpi
    USE cell_base,          ONLY : bg, at, alat, tpiba, tpiba2, omega
    USE gvect,              ONLY : ngm, g
    USE gvecw,              ONLY : gcutw
    USE mp_global,          ONLY : intra_pool_comm
    USE mp,                 ONLY : mp_sum
    USE io_global,          ONLY : stdout
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: skip_eps
    !
    INTEGER, PARAMETER :: nrad = 4000      ! Simpson points for the radial integral (even)
    INTEGER, PARAMETER :: nth  = 360       ! angular points (midpoint)
    INTEGER  :: iq1, iq2, ig, ir, ith
    REAL(DP) :: c, lz, area, beta, alpha_2d(2,2), r0_mat(2,2)
    REAL(DP) :: xq(3), q(2), qn, qq, f, r0q, sum_b, sum_s
    REAL(DP) :: int_b, int_s, qmax, dq, th, u(2), r0u, rad, w, pref, sheet
    REAL(DP) :: exxdiv_3d, exxdiv_eps_3d
    !
    c    = at(3,3)*alat
    lz   = 0.5_DP*c
    area = omega/c
    beta = 10._DP/(gcutw*tpiba2)
    !
    alpha_2d(:,:) = 0._DP
    IF (.NOT. skip_eps) THEN
       alpha_2d(:,:) = eps_mat(1:2,1:2)
       alpha_2d(1,1) = alpha_2d(1,1) - 1._DP
       alpha_2d(2,2) = alpha_2d(2,2) - 1._DP
       alpha_2d(:,:) = alpha_2d(:,:) * c / fpi
    ENDIF
    r0_mat(:,:) = tpi*alpha_2d(:,:)
    !
    ! ... grid sum of the auxiliary functions over the q-mesh and the in-plane G
    sum_b = 0._DP
    sum_s = 0._DP
    DO iq1 = 0, nq1-1
      DO iq2 = 0, nq2-1
        xq(:) = bg(:,1)*DBLE(iq1)/DBLE(nq1) + bg(:,2)*DBLE(iq2)/DBLE(nq2)
        DO ig = 1, ngm
          IF ( ABS(g(3,ig)) > 1.D-8 ) CYCLE                ! in-plane G only
          q(1:2) = ( xq(1:2) + g(1:2,ig) ) * tpiba          ! bohr^-1
          qq = q(1)**2 + q(2)**2
          IF ( qq < 1.D-12 ) CYCLE                          ! q+G = 0
          qn  = SQRT(qq)
          f   = fpi*e2*lz*EXP(-beta*qq)/qn
          r0q = DOT_PRODUCT( q, MATMUL(r0_mat, q) )/qn      ! r0(qhat) * q
          sum_b = sum_b + f
          sum_s = sum_s + f/(1._DP + r0q)
        ENDDO
      ENDDO
    ENDDO
    CALL mp_sum( sum_b, intra_pool_comm )
    CALL mp_sum( sum_s, intra_pool_comm )
    !
    ! ... integrals of the auxiliary functions over the plane
    int_b = fpi*e2*lz * tpi * 0.5_DP*SQRT(pi/beta)
    qmax  = SQRT(60._DP/beta)
    dq    = qmax/DBLE(nrad)
    int_s = 0._DP
    DO ith = 1, nth
      th  = tpi*(DBLE(ith)-0.5_DP)/DBLE(nth)
      u   = (/ COS(th), SIN(th) /)
      r0u = DOT_PRODUCT( u, MATMUL(r0_mat, u) )
      rad = 0._DP
      DO ir = 0, nrad
        w = 2._DP
        IF (MOD(ir,2) == 1) w = 4._DP
        IF (ir == 0 .OR. ir == nrad) w = 1._DP
        qn  = DBLE(ir)*dq
        rad = rad + w*EXP(-beta*qn*qn)/(1._DP + r0u*qn)
      ENDDO
      int_s = int_s + rad*dq/3._DP
    ENDDO
    int_s = fpi*e2*lz * int_s*tpi/DBLE(nth)
    !
    pref  = DBLE(nq1*nq2)*area/tpi**2
    sheet = tpi*e2*lz**2
    !
    exxdiv_3d     = exxdiv
    exxdiv_eps_3d = exxdiv_eps
    exxdiv        = - ( pref*int_b - sum_b - sheet )
    exxdiv_eps    = - ( pref*int_s - sum_s - sheet )
    !
    WRITE(stdout,'(/,5X, "INFO: 2D cutoff q+G=0 treatment (charged sheet + 2D Gygi-Baldereschi)")')
    WRITE(stdout,'(  5X, "INFO: lz = c/2 [bohr]     ", 3x, 1F12.6)') lz
    WRITE(stdout,'(  5X, "INFO: GB width beta [bohr^2]", 1x, 1F12.6)') beta
    WRITE(stdout,'(  5X, "INFO: alpha_2D [bohr]     ", 3x, 2F12.6)') alpha_2d(:,1)
    WRITE(stdout,'(  5X, "                          ", 3x, 2F12.6)') alpha_2d(:,2)
    WRITE(stdout,'(  5X, "INFO: r0 = 2pi alpha_2D   ", 3x, 1F12.6, " bohr (in-plane average)")') &
                 0.5_DP*(r0_mat(1,1)+r0_mat(2,2))
    WRITE(stdout,'(  5X, "INFO: F0 bare, screened   ", 3x, 2F20.10)') pref*int_b-sum_b, pref*int_s-sum_s
    WRITE(stdout,'(  5X, "INFO: 3D GB estimates     ", 3x, 2ES15.5, "  (NOT used)")') exxdiv_3d, exxdiv_eps_3d
    IF (skip_eps) THEN
       ! only reachable for calculations that do not use exxdiv_eps (setup_coulomb stops
       ! for screen and cc without eps.dat)
       WRITE(stdout,'(5X, "INFO: no dielectric information: alpha_2D = 0 (screened q+G=0 term", &
                    & " not used in this calculation)")')
    ELSE IF (ABS(eps_mat(3,3)-eps_mat(1,1)) < 1.D-8 .AND. ABS(eps_mat(1,2)) < 1.D-8 .AND. &
             ABS(eps_mat(1,1)-eps_mat(2,2)) < 1.D-8) THEN
       WRITE(stdout,'(5X, "WARNING: isotropic eps used for a 2D system: provide the supercell", &
                    & " tensor in eps.dat (only the in-plane block is used)")')
    ENDIF
    !
  END SUBROUTINE divergence_2d
  !
  !---------------------------------------------------------------------
  SUBROUTINE read_eps (eps_mat) 
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL (DP) :: eps_mat (3,3)
    INTEGER   :: i, j, ios
    !
    OPEN (765, file = 'eps.dat')
    DO j = 1, 3
      READ (765, *, IOSTAT=ios) (eps_mat(j,i), i =1,3)
      IF (ios /= 0) THEN 
       WRITE(stdout,'(/, 5X,A, I5)') "ERROR: Somethng wrong reading eps.dat", ios
       STOP
      ENDIF
    ENDDO
    CLOSE (765)
    !
    RETURN
    !
  END SUBROUTINE read_eps
  !
END MODULE
