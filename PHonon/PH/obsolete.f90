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
