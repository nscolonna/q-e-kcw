!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  This file contains a few obsolete routines that are still used
!  (but shouldn't) in some parts of the phonon code.
!
!-----------------------------------------------------------------------
subroutine smallgq (xq, at, bg, s, nsym, irgq, nsymq, irotmq, &
     minus_q, gi, gimq)
  !-----------------------------------------------------------------------
  !! -- OBSOLETE (kept for compatibility) --  
  !! This routine selects, among the symmetry matrices of the point group
  !! of a crystal, the symmetry operations which leave q unchanged.
  !! Furthermore it checks if one of the matrices send q <-> -q+G. In
  !! this case \(\text{minus_q}\) is set TRUE.
  !
  !! Revised   2 Sept. 1995 by Andrea Dal Corso
  !! Modified 22 April 1997 by SdG: \(\text{minus_q}\) is sought also among
  !!   sym.op. such that Sq=q+G (i.e. the case \(q=-q+G\) is dealt with).
  !
  !
  !  The dummy variables
  !
  USE kinds, only : DP
  implicit none

  real(DP), parameter :: accep=1.e-5_dp
  real(DP) :: bg (3, 3), at (3, 3), xq (3), gi (3, 48), gimq (3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the q point of the crystal
  ! output: the G associated to a symmetry:[S(irotq)*q - q]
  ! output: the G associated to:  [S(irotmq)*q + q]

  integer :: s (3, 3, 48), irgq (48), irotmq, nsymq, nsym
  ! input: the symmetry matrices
  ! output: the symmetry of the small group
  ! output: op. symmetry: s_irotmq(q)=-q+G
  ! output: dimension of the small group of q
  ! input: dimension of the point group

  logical :: minus_q
  ! input: .t. if sym.ops. such that Sq=-q+G are searched for
  ! output: .t. if such a symmetry has been found

  real(DP) :: wrk (3), aq (3), raq (3), zero (3)
  ! additional space to compute gi and gimq
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: isym, ipol, jpol
  ! counter on symmetry operations
  ! counter on polarizations
  ! counter on polarizations

  logical :: look_for_minus_q, eqvect
  ! .t. if sym.ops. such that Sq=-q+G are searched for
  ! logical function, check if two vectors are equal
  !
  !  Set to zero some variables and transform xq to the crystal basis
  !
  look_for_minus_q = minus_q
  !
  minus_q = .false.
  zero = 0.d0
  gi   = 0.d0
  gimq = 0.d0
  aq = xq
  call cryst_to_cart (1, aq, at, - 1)
  !
  !   test all symmetries to see if the operation S sends q in q+G ...
  !
  nsymq = 0
  do isym = 1, nsym
     raq = 0.d0
     do ipol = 1, 3
        do jpol = 1, 3
           raq (ipol) = raq (ipol) + DBLE (s (ipol, jpol, isym) ) * &
                aq (jpol)
        enddo
     enddo
     if (eqvect (raq, aq, zero, accep) ) then
        nsymq = nsymq + 1
        irgq (nsymq) = isym
        do ipol = 1, 3
           wrk (ipol) = raq (ipol) - aq (ipol)
        enddo
        call cryst_to_cart (1, wrk, bg, 1)
        gi (:, nsymq) = wrk (:)
        !
        !   ... and in -q+G
        !
        if (look_for_minus_q.and..not.minus_q) then
           raq (:) = - raq(:)
           if (eqvect (raq, aq, zero, accep) ) then
              minus_q = .true.
              irotmq = isym
              do ipol = 1, 3
                 wrk (ipol) = - raq (ipol) + aq (ipol)
              enddo
              call cryst_to_cart (1, wrk, bg, 1)
              gimq (:) = wrk (:)
           endif
        endif
     endif
  enddo
  !
  ! if xq=(0,0,0) minus_q always apply with the identity operation
  !
  if (xq (1) == 0.d0 .and. xq (2) == 0.d0 .and. xq (3) == 0.d0) then
     minus_q = .true.
     irotmq = 1
     gimq = 0.d0
  endif
  !
  return
end subroutine smallgq

!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  !-----------------------------------------------------------------------
  !! -- OBSOLETE (kept for compatibility) --  
  !! This routine computes the vector \(\text{rtau}\) which contains for each
  !! atom and each rotation the vector \(S\tau_a - \tau_b\), where
  !! \(b\) is the rotated \(a\) atom, given by the array \(\text{irt}. These 
  !! \(\text{rtau}\) are non zero only if fractional translations are 
  !! present.
  !
  USE kinds, ONLY : DP
  implicit none
  !
  !     first the dummy variables
  !
  integer, intent(in) :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! nsym: number of symmetries of the point group
  ! s:    matrices of symmetry operations
  ! nat : number of atoms in the unit cell
  ! irt(n,m) = transformed of atom m for symmetry n
  real(DP), intent(in) :: at (3, 3), bg (3, 3), tau (3, nat)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  ! tau: coordinates of the atoms
  logical, intent(in) :: sym (nsym)
  ! sym(n)=.true. if operation n is a symmetry
  real(DP), intent(out):: rtau (3, 48, nat)
  ! rtau: the direct translations
  !
  !    here the local variables
  !
  integer :: na, nb, isym, ipol
  ! counters on: atoms, symmetry operations, polarization
  real(DP) , allocatable :: xau (:,:)
  real(DP) :: ft (3)
  !
  allocate (xau(3,nat))
  !
  !   compute the atomic coordinates in crystal axis, xau
  !
  do na = 1, nat
     do ipol = 1, 3
        xau (ipol, na) = bg (1, ipol) * tau (1, na) + &
                         bg (2, ipol) * tau (2, na) + &
                         bg (3, ipol) * tau (3, na)
     enddo
  enddo
  !
  !    for each symmetry operation, compute the atomic coordinates
  !    of the rotated atom, ft, and calculate rtau = Stau'-tau
  !
  rtau(:,:,:) = 0.0_dp
  do isym = 1, nsym
     if (sym (isym) ) then
        do na = 1, nat
           nb = irt (isym, na)
           do ipol = 1, 3
              ft (ipol) = s (1, ipol, isym) * xau (1, na) + &
                          s (2, ipol, isym) * xau (2, na) + &
                          s (3, ipol, isym) * xau (3, na) - xau (ipol, nb)
           enddo
           do ipol = 1, 3
              rtau (ipol, isym, na) = at (ipol, 1) * ft (1) + &
                                      at (ipol, 2) * ft (2) + &
                                      at (ipol, 3) * ft (3)
           enddo
        enddo
     endif
  enddo
  !
  !    deallocate workspace
  !
  deallocate(xau)
  return
end subroutine sgam_ph
!

!
! Copyright (C) 2006-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_mode_sym (u, w2, at, bg, tau, nat, nsym, sr, irt, xq, &
     rtau, amass, ntyp, ityp, flag, lri, lmolecule, nspin_mag,        &
     name_rap_mode, num_rap_mode)
  !
  !! -- OBSOLETE (kept for compatibility) --  
  !! This subroutine finds the irreducible representations which give
  !! the transformation properties of eigenvectors of the dynamical
  !! matrix. It does NOT work at zone border in non symmorphic space groups.
  !! if flag=1 the true displacements are given in input, otherwise the
  !! eigenvalues of the dynamical matrix are given.
  !
  !
  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, RY_TO_CMM1
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, which_irr, &
       char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_is, ONLY : gname_is
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::             &
       nat,         &
       nsym,        &
       flag,        &
       ntyp,        &
       nspin_mag,   &
       ityp(nat),   &
       irt(48,nat)

  REAL(DP), INTENT(IN) ::   &
       at(3,3),        &
       bg(3,3),        &
       xq(3),          &
       tau(3,nat),     &
       rtau(3,48,nat), &
       amass(ntyp),    &
       w2(3*nat),      &
       sr(3,3,48)

  COMPLEX(DP), INTENT(IN) ::  &
       u(3*nat, 3*nat)       ! The eigenvectors or the displacement pattern
  LOGICAL, INTENT(IN) :: lri      ! if .true. print the Infrared/Raman flag
  LOGICAL, INTENT(IN) :: lmolecule ! if .true. these are eigenvalues of an
                                   ! isolated system

  CHARACTER(15), INTENT(OUT) :: name_rap_mode( 3 * nat )
  INTEGER, INTENT(OUT) :: num_rap_mode ( 3 * nat )

  REAL(DP), PARAMETER :: eps=1.d-5

  INTEGER ::      &
       ngroup, &   ! number of different frequencies groups
       nmodes, &   ! number of modes
       imode, imode1, igroup, dim_rap, nu_i,  &
       irot, irap, iclass, mu, na, i

  INTEGER, ALLOCATABLE :: istart(:)

  COMPLEX(DP) :: times              ! safe dimension
  ! in case of accidental degeneracy
  REAL(DP), ALLOCATABLE :: w1(:)
  COMPLEX(DP), ALLOCATABLE ::  rmode(:,:), trace(:,:), z(:,:)
  LOGICAL :: is_linear
  CHARACTER(3) :: cdum
  INTEGER :: counter, counter_s
  !
  !    Divide the modes on the basis of the mode degeneracy.
  !
  nmodes=3*nat

  ALLOCATE(istart(nmodes+1))
  ALLOCATE(z(nmodes,nmodes))
  ALLOCATE(w1(nmodes))
  ALLOCATE(rmode(nmodes,nmodes))
  ALLOCATE(trace(48,nmodes))

  IF (flag==1) THEN
     !
     !  Find the eigenvalues of the dynmaical matrix
     !  Note that amass is in amu; amu_ry converts it to Ry au
     !
     DO nu_i = 1, nmodes
        DO mu = 1, nmodes
           na = (mu - 1) / 3 + 1
           z (mu, nu_i) = u (mu, nu_i) * SQRT (amu_ry*amass (ityp (na) ) )
        END DO
     END DO
  ELSE
     z=u
  ENDIF

  DO imode=1,nmodes
     w1(imode)=SIGN(SQRT(ABS(w2(imode)))*RY_TO_CMM1,w2(imode))
  ENDDO

  ngroup=1
  istart(ngroup)=1
  imode1=1
  IF (lmolecule) THEN
     istart(ngroup)=7
     imode1=6
     IF(is_linear(nat,tau)) istart(ngroup)=6
  ENDIF
  DO imode=imode1+1,nmodes
     IF (ABS(w1(imode)-w1(imode-1)) > 5.0d-2) THEN
        ngroup=ngroup+1
        istart(ngroup)=imode
     END IF
  END DO
  istart(ngroup+1)=nmodes+1
  !
  !  Find the character of one symmetry operation per class
  !
  DO igroup=1,ngroup
     dim_rap=istart(igroup+1)-istart(igroup)
     DO iclass=1,nclass
        irot=elem(1,iclass)
        trace(iclass,igroup)=(0.d0,0.d0)
        DO i=1,dim_rap
           nu_i=istart(igroup)+i-1
           CALL rotate_mod(z,rmode,sr(1,1,irot),irt,rtau,xq,nat,irot)
           trace(iclass,igroup)=trace(iclass,igroup) + &
                   dot_product(z(1:3*nat,nu_i),rmode(1:3*nat,nu_i))
        END DO
!              write(6,*) igroup,iclass, trace(iclass,igroup)
     END DO
  END DO
  !
  !  And now use the character table to identify the symmetry representation
  !  of each group of modes
  !
  IF (nspin_mag==4) THEN
     IF (flag==1) WRITE(stdout,  &
          '(/,5x,"Mode symmetry, ",a11," [",a11,"] magnetic point group:",/)') &
          gname, gname_is
  ELSE
     IF (flag==1) WRITE(stdout,'(/,5x,"Mode symmetry, ",a11," point group:",/)') gname
  END IF
  num_rap_mode=-1
  counter=1
  DO igroup=1,ngroup
     IF (ABS(w1(istart(igroup)))<1.d-3) CYCLE
     DO irap=1,nclass
        times=(0.d0,0.d0)
        DO iclass=1,nclass
           times=times+CONJG(trace(iclass,igroup))*char_mat(irap, &
                which_irr(iclass))*nelem(iclass)
           !         write(6,*) igroup, irap, iclass, which_irr(iclass)
        ENDDO
        times=times/nsym
        cdum="   "
        IF (lri) cdum=ir_ram(irap)
        IF ((ABS(NINT(DBLE(times))-DBLE(times)) > 1.d-4).OR. &
             (ABS(AIMAG(times)) > eps) ) THEN
           IF (flag==1) WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x, "-->   ?")') &
                istart(igroup), istart(igroup+1)-1, w1(istart(igroup))
        ENDIF

        IF (ABS(times) > eps) THEN
           IF (ABS(NINT(DBLE(times))-1.d0) < 1.d-4) THEN
              IF (flag==1) WRITE(stdout,'(5x, "omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",a19)') &
                   istart(igroup), istart(igroup+1)-1, w1(istart(igroup)), &
                   name_rap(irap)//" "//cdum
              name_rap_mode(igroup)=name_rap(irap)
              counter_s=counter
              DO imode=counter_s, counter_s+NINT(DBLE(char_mat(irap,1)))-1
                 IF (imode <= 3*nat) num_rap_mode(imode) = irap
                 counter=counter+1
              ENDDO
           ELSE
              IF (flag==1) WRITE(stdout,'(5x,"omega(",i3," -",i3,") = ",f12.1,2x,"[cm-1]",3x,"--> ",i3,a19)') &
                   istart(igroup), istart(igroup+1)-1, &
                   w1(istart(igroup)), NINT(DBLE(times)), &
                   name_rap(irap)//" "//cdum
              name_rap_mode(igroup)=name_rap(irap)
              counter_s=counter
              DO imode=counter_s, counter_s+NINT(DBLE(times))*&
                                            NINT(DBLE(char_mat(irap,1)))-1
                 IF (imode <= 3 * nat) num_rap_mode(imode) = irap
                 counter=counter+1
              ENDDO
           END IF
        END IF
     END DO
  END DO
  IF (flag==1) WRITE( stdout, '(/,1x,74("*"))')

  DEALLOCATE(trace)
  DEALLOCATE(z)
  DEALLOCATE(w1)
  DEALLOCATE(rmode)
  DEALLOCATE(istart)

  RETURN
END SUBROUTINE find_mode_sym
