!
! Copyright (C) 2001-2007 Quantum ESPRESSO PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine add_for_charges (ik, uact)
  !----------===============-----------------------------------------------
  !! This subroutine calculates: \(\frac{dS}{du} P_c [x, H-eS] |\psi\rangle\)
  !

  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE cell_base, ONLY : tpiba
  USE gvect, ONLY : g
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE klist, ONLY : xk, ngk, igk_k
  USE uspp, ONLY : nkb, qq_nt, qq_so, vkb, ofsbeta
  USE wvfct, ONLY : npwx, nbnd
  USE becmod, ONLY: calbec, bec_type, allocate_bec_type, deallocate_bec_type
  USE noncollin_module, ONLY : noncolin, npol, lspinorb
  USE uspp_param, only: nh, nhm
  USE eqv, ONLY : dvpsi, dpsi
  USE control_lr, ONLY : lgamma
  USE qpoint,  ONLY : ikks
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik
  !! input: the k point
  integer :: mode
  ! input: the actual perturbation
  complex(DP) :: uact (3 * nat)
  !! input: the pattern of displacements
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol, is, js, ijs, npw
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! the point k+q
  ! counter on G vectors
  ! auxiliary counter on G vectors
  ! counter on atomic types
  ! counter on bands
  ! auxiliary variable for counting
  ! counter on becp functions
  ! counter on becp functions
  ! counter on n index
  ! counter on m index
  ! counter on polarizations

  real(DP), parameter :: eps = 1.d-12

  complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:)
  complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  ! temporary arrays for optimization
  complex(DP) :: temp_ps1, temp_ps2(3)
  complex(DP) :: temp_ps1_nc(npol), temp_ps2_nc(npol,3)
  complex(DP) :: temp_ps2_nc_1, temp_ps2_nc_2, temp_ps2_k
  ! small buffers for loop optimization (allocated once with max size)
  complex(DP) :: alphapp_buf_nc(nhm, npol), bedp_buf_nc(nhm, npol)
  complex(DP) :: alphapp_buf_k(nhm), bedp_buf_k(nhm)
  complex(DP) :: qq_so_buf(nhm, npol*npol), qq_nt_buf(nhm)
  ! the scalar product
  ! the scalar product
  ! a mesh space for psi
  TYPE(bec_type) :: bedp, alphapp(3)
  complex(DP), allocatable :: aux1(:,:)

  logical :: ok
  ! used to save time

  allocate (aux ( npwx))
  allocate (aux1( npwx*npol, nbnd))
  CALL allocate_bec_type(nkb,nbnd,bedp)
  DO ipol=1,3
     CALL allocate_bec_type(nkb,nbnd,alphapp(ipol))
  ENDDO
  IF (noncolin) THEN
     allocate (ps1_nc ( nkb, npol, nbnd))
     allocate (ps2_nc ( nkb, npol, nbnd , 3))
  ELSE
     allocate (ps1 ( nkb , nbnd))
     allocate (ps2 ( nkb , nbnd , 3))
  ENDIF
  if (lgamma) then
     ikk = ikks(ik)
     ikq = ikk
     npw =ngk(ikk)
  else
     call infomsg ('add_for_charges', 'called for lgamma .eq. false')
  endif
  if (lsda) current_spin = isk (ikk)
  !
  !   we first compute the coefficients of the vectors
  !
  if (noncolin) then
     ps1_nc   = (0.d0, 0.d0)
     ps2_nc   = (0.d0, 0.d0)
     bedp%nc = (0.d0,0.d0)
     DO ipol=1,3
        alphapp(ipol)%nc = (0.d0,0.d0)
     END DO
  else
     ps1   = (0.d0, 0.d0)
     ps2   = (0.d0, 0.d0)
     bedp%k = (0.d0,0.d0)
     DO ipol=1,3
        alphapp(ipol)%k = (0.d0,0.d0)
     END DO
  endif
  aux1  = (0.d0, 0.d0)

  !
  ! first we calculate the products of the beta functions with dpsi
  !
  CALL calbec (npw, vkb, dpsi, bedp)
  do ipol = 1, 3
     aux1=(0.d0,0.d0)
     do ibnd = 1, nbnd
        do ig = 1, npw
           aux1 (ig, ibnd) = dpsi(ig,ibnd) *           &
                tpiba * (0.d0,1.d0) *                  &
                ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
        enddo
        if (noncolin) then
           do ig = 1, npw
              aux1 (ig+npwx, ibnd) = dpsi(ig+npwx,ibnd) *           &
                   tpiba * (0.d0,1.d0) *                  &
                  ( xk(ipol,ikk) + g(ipol,igk_k(ig,ikk)) )
           enddo
        endif
     enddo
     CALL calbec ( npw, vkb, aux1, alphapp(ipol) )
  enddo


  do na = 1, nat
     nt = ityp (na)
     ijkb0=ofsbeta(na)
     mu = 3 * (na - 1)
     if ( abs (uact (mu + 1) ) + &
          abs (uact (mu + 2) ) + &
          abs (uact (mu + 3) ) > eps) then
        do ih = 1, nh (nt)
           ikb = ijkb0 + ih
           do ipol = 1, 3
              do ibnd = 1, nbnd
                 ! Populate small buffers for this iteration
                 qq_nt_buf(1:nh(nt)) = qq_nt(ih, 1:nh(nt), nt)
                 
                 if (noncolin) then
                    do jh = 1, nh(nt)
                       jkb = ijkb0 + jh
                       alphapp_buf_nc(jh, 1:npol) = alphapp(ipol)%nc(jkb, 1:npol, ibnd)
                       bedp_buf_nc(jh, 1:npol) = bedp%nc(jkb, 1:npol, ibnd)
                    enddo
                    if (lspinorb) then
                       do jh = 1, nh(nt)
                          qq_so_buf(jh, 1:npol*npol) = qq_so(ih, jh, 1:npol*npol, nt)
                       enddo
                    endif
                 else
                    do jh = 1, nh(nt)
                       jkb = ijkb0 + jh
                       alphapp_buf_k(jh) = alphapp(ipol)%k(jkb, ibnd)
                       bedp_buf_k(jh) = bedp%k(jkb, ibnd)
                    enddo
                 endif
                 
                 ! Initialize temp arrays
                 if (noncolin) then
                    temp_ps1_nc = (0.d0, 0.d0)
                    temp_ps2_nc = (0.d0, 0.d0)
                 else
                    temp_ps1 = (0.d0, 0.d0)
                    temp_ps2 = (0.d0, 0.d0)
                 endif
                 
                 do jh = 1, nh (nt)
                    if (noncolin) then
                       if (lspinorb) then
                          ijs=0
                          DO is=1,npol
                             DO js=1,npol
                                ijs=ijs+1
                                temp_ps1_nc(is) = temp_ps1_nc(is) + &
                                (qq_so_buf(jh,ijs) *              &
                                alphapp_buf_nc(jh,js))*         &
                                uact (mu + ipol)
                                temp_ps2_nc(is,ipol) = temp_ps2_nc(is,ipol) + &
                                (qq_so_buf(jh,ijs) *              &
                                 bedp_buf_nc(jh,js))*(0.d0,-1.d0)* &
                                 uact (mu + ipol) * tpiba
                             ENDDO
                          ENDDO
                       else
                          do is=1,npol
                             temp_ps1_nc(is) = temp_ps1_nc(is) + &
                                 qq_nt_buf(jh) *                     &
                                 alphapp_buf_nc(jh,is) *     &
                                 uact (mu + ipol)
                             temp_ps2_nc(is,ipol) = temp_ps2_nc(is,ipol) + &
                                 qq_nt_buf(jh) * (0.d0, -1.d0) *     &
                                 bedp_buf_nc(jh,is) *             &
                                 uact (mu + ipol) * tpiba
                          end do
                       endif
                    else
                       temp_ps1 = temp_ps1 + qq_nt_buf(jh)*alphapp_buf_k(jh)* &
                            uact (mu + ipol)
                       temp_ps2(ipol) = temp_ps2(ipol) + qq_nt_buf(jh) * (0.d0, -1.d0) * &
                             bedp_buf_k(jh) *uact (mu + ipol) * tpiba
                    endif
                 enddo
                 
                 ! Assign temp values to main arrays
                 if (noncolin) then
                    do is=1,npol
                       ps1_nc(ikb,is,ibnd) = ps1_nc(ikb,is,ibnd) + temp_ps1_nc(is)
                       ps2_nc(ikb,is,ibnd,ipol) = ps2_nc(ikb,is,ibnd,ipol) + temp_ps2_nc(is,ipol)
                    enddo
                 else
                    ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + temp_ps1
                    ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + temp_ps2(ipol)
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  if (nkb.gt.0) then
     if (noncolin) then
        call zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
     else
        call zgemm ('N', 'N', npw, nbnd*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
!        dvpsi = matmul(vkb,ps1) + dvpsi
     endif
  endif
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        if (noncolin) then
           ok = ANY(ABS(ps2_nc(ikb, 1, 1:nbnd, ipol)) > eps) .OR. &
                ANY(ABS(ps2_nc(ikb, 2, 1:nbnd, ipol)) > eps)
        else
           ok = ANY(ABS(ps2(ikb, 1:nbnd, ipol)) > eps)
        endif
        if (ok) then
           do ig = 1, npw
              igg = igk_k (ig,ikq)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, ikq) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              if (noncolin) then
                 temp_ps2_nc_1 = ps2_nc(ikb, 1, ibnd, ipol)
                 temp_ps2_nc_2 = ps2_nc(ikb, 2, ibnd, ipol)
                 do ig = 1, npw
                    dvpsi(ig,      ibnd) = temp_ps2_nc_1*aux(ig) + dvpsi(ig,      ibnd)
                    dvpsi(ig+npwx, ibnd) = temp_ps2_nc_2*aux(ig) + dvpsi(ig+npwx, ibnd)
                 enddo
              else
                 temp_ps2_k = ps2(ikb, ibnd, ipol)
                 do ig = 1, npw
                    dvpsi(ig, ibnd) = temp_ps2_k*aux(ig) + dvpsi(ig, ibnd)
                 enddo
              endif
           enddo
        endif
     enddo
  enddo
!
!    Now dvpsi contains dS/du x |psi>
!

  deallocate (aux)
  deallocate (aux1)
  IF (noncolin) THEN
     deallocate (ps1_nc)
     deallocate (ps2_nc)
  ELSE
     deallocate (ps1)
     deallocate (ps2)
  END IF
  CALL deallocate_bec_type(bedp)
  DO ipol=1,3
     CALL deallocate_bec_type(alphapp(ipol))
  END DO

  return
end subroutine add_for_charges
