!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us_only (ik, uact, becp1, alphap)
   !----------------------------------------------------------------------
   !! This routine calculates \(\text{dV_bare}/\text{dtau}\cdot\text{psi}\)
   !! for one perturbation with a given q.  
   !! The displacements are described by a vector uact.  
   !! The result is stored in \(\text{dvpsi}\). The routine is called for
   !! each k-point and for each pattern u. It computes simultaneously all
   !! the bands.  
   !! This routine implements Eq. (B29) of PRB 64, 235118 (2001).
   !! Only the contribution of the nonlocal potential is calculated here;
   !! both norm-conserving term and ultrasoft correction are calculated here.
   !
   !
   USE kinds, only : DP
   USE cell_base, ONLY : tpiba
   USE gvect,     ONLY : g
   USE klist,     ONLY : xk, ngk, igk_k
   USE ions_base, ONLY : nat, ityp, ntyp => nsp
   USE lsda_mod,  ONLY : lsda, current_spin, isk, nspin
   USE wvfct,     ONLY : nbnd, npwx, et
   USE noncollin_module, ONLY : noncolin, npol, lspinorb
   USE uspp, ONLY: okvan, nkb, vkb, ofsbeta
   USE uspp_param, ONLY: nh, nhm
   USE phus,      ONLY : int1, int1_nc, int2, int2_so

   USE qpoint,     ONLY : nksq, ikks, ikqs
   USE becmod,     ONLY : bec_type
   USE eqv,        ONLY : dvpsi
   USE control_lr, ONLY : lgamma

   implicit none
   !
   integer :: ik
   !! input: the k point
   complex(DP) :: uact(3*nat)
   !! input: the pattern of displacements
   TYPE(bec_type) :: becp1(nksq), alphap(3,nksq)
   !
   !   ... local variables
   !

   integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
        ikb, jkb, ih, jh, ipol, is, js, ijs, npwq
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

   complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:), deff_nc(:,:,:,:)
   real(DP), allocatable :: deff(:,:,:)
   complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
   ! work space
   complex(DP), allocatable :: temp_ps1 (:), temp_ps2 (:,:)
   complex(DP), allocatable :: temp_ps1_nc (:,:), temp_ps2_nc (:,:,:)
   ! temporary buffers

   logical :: ok

   call start_clock ('dvqpsi_us_on')
   if (noncolin) then
      allocate (ps1_nc(nkb , npol, nbnd))
      allocate (ps2_nc(nkb , npol, nbnd , 3))
      allocate (deff_nc(nhm, nhm, nat, nspin))
      allocate (temp_ps1_nc(nkb, npol))
      allocate (temp_ps2_nc(nkb, npol, 3))
   else
      allocate (ps1 ( nkb , nbnd))
      allocate (ps2 ( nkb , nbnd , 3))
      allocate (deff(nhm, nhm, nat))
      allocate (temp_ps1(nkb))
      allocate (temp_ps2(nkb, 3))
   end if
   allocate (aux ( npwx))
   ikk = ikks(ik)
   ikq = ikqs(ik)
   if (lsda) current_spin = isk (ikk)
   !
   !   we first compute the coefficients of the vectors
   !
   if (noncolin) then
      ps1_nc(:,:,:)   = (0.d0, 0.d0)
      ps2_nc(:,:,:,:) = (0.d0, 0.d0)
   else
      ps1(:,:)   = (0.d0, 0.d0)
      ps2(:,:,:) = (0.d0, 0.d0)
   end if
   do ibnd = 1, nbnd
      ! Initialize temporary buffers
      if (noncolin) then
         temp_ps1_nc(:,:) = (0.d0, 0.d0)
         temp_ps2_nc(:,:,:) = (0.d0, 0.d0)
      else
         temp_ps1(:) = (0.d0, 0.d0)
         temp_ps2(:,:) = (0.d0, 0.d0)
      end if
      
      IF (noncolin) THEN
         CALL compute_deff_nc(deff_nc,et(ibnd,ikk))
      ELSE
         CALL compute_deff(deff,et(ibnd,ikk))
      ENDIF
      do na = 1, nat
         ijkb0 = ofsbeta(na) 
         nt = ityp(na) 
         mu = 3 * (na - 1)
         ! First loop: deff calculation
         if ( abs (uact (mu + 1) ) + &
              abs (uact (mu + 2) ) + &
              abs (uact (mu + 3) ) > eps) then
            do ih = 1, nh (nt)
               ikb = ijkb0 + ih
               do jh = 1, nh (nt)
                  jkb = ijkb0 + jh
                  do ipol = 1, 3
                     IF (noncolin) THEN
                        ijs=0
                        DO is=1,npol
                           DO js=1,npol
                              ijs=ijs+1
                              temp_ps1_nc(ikb,is)=temp_ps1_nc(ikb,is) +  &
                                 deff_nc(ih,jh,na,ijs) * &
                                 alphap(ipol, ik)%nc(jkb,js,ibnd)* &
                                  uact(mu + ipol)
                              temp_ps2_nc(ikb,is,ipol)=               &
                                     temp_ps2_nc(ikb,is,ipol)+        &
                                     deff_nc(ih,jh,na,ijs) *          &
                                     becp1(ik)%nc(jkb,js,ibnd) *      &
                                     (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                           END DO
                        END DO
                     ELSE
                        temp_ps1 (ikb) = temp_ps1 (ikb) +      &
                                   deff(ih, jh, na) *            &
                           alphap(ipol, ik)%k(jkb, ibnd) * uact (mu + ipol)
                        temp_ps2 (ikb, ipol) = temp_ps2 (ikb, ipol) +&
                             deff(ih,jh,na)*becp1(ik)%k (jkb, ibnd) *  &
                             (0.0_DP,-1.0_DP) * uact (mu + ipol) * tpiba
                     ENDIF
                     IF (okvan) THEN
                        IF (noncolin) THEN
                           ijs=0
                           DO is=1,npol
                              DO js=1,npol
                                 ijs=ijs+1
                                 temp_ps1_nc(ikb,is)=temp_ps1_nc(ikb,is)+ &
                                    int1_nc(ih,jh,ipol,na,ijs) *     &
                                    becp1(ik)%nc(jkb,js,ibnd)*uact(mu+ipol)
                              END DO
                           END DO
                        ELSE
                           temp_ps1 (ikb) = temp_ps1 (ikb) + &
                             (int1 (ih, jh, ipol,na, current_spin) * &
                             becp1(ik)%k (jkb, ibnd) ) * uact (mu +ipol)
                        END IF
                     END IF
                  enddo ! ipol
               enddo ! jh
            enddo ! ih
         END IF  ! uact>0
      enddo  ! na
      ! Second loop: okvan nb loop
      if (okvan) then
         do nb = 1, nat
            nu = 3 * (nb - 1)
            if ( abs (uact (nu + 1) ) + &
                 abs (uact (nu + 2) ) + &
                 abs (uact (nu + 3) ) > eps) then
               do na = 1, nat
                  ijkb0 = ofsbeta(na) 
                  nt = ityp(na) 
                  do ih = 1, nh (nt)
                     ikb = ijkb0 + ih
                     do jh = 1, nh (nt)
                        jkb = ijkb0 + jh
                        do ipol = 1, 3
                           IF (noncolin) THEN
                              IF (lspinorb) THEN
                                 ijs=0
                                 DO is=1,npol
                                    DO js=1,npol
                                       ijs=ijs+1
                                       temp_ps1_nc(ikb,is)= &
                                                 temp_ps1_nc(ikb,is)+ &
                                       int2_so(ih,jh,ipol,nb,na,ijs)* &
                                        becp1(ik)%nc(jkb,js,ibnd)*uact(nu+ipol)
                                    END DO
                                 END DO
                              ELSE
                                 DO is=1,npol
                                    temp_ps1_nc(ikb,is)=temp_ps1_nc(ikb,is)+ &
                                       int2(ih,jh,ipol,nb,na) * &
                                       becp1(ik)%nc(jkb,js,ibnd)*uact(nu+ipol)
                                 END DO
                              END IF
                           ELSE
                              temp_ps1 (ikb) = temp_ps1 (ikb) + &
                                  (int2 (ih, jh, ipol, nb, na) * &
                                   becp1(ik)%k (jkb, ibnd) ) * uact (nu + ipol)
                           END IF
                        enddo ! ipol
                     enddo ! jh
                  enddo ! ih
               enddo  ! na
            END IF  ! uact>0
         enddo  ! nb
      endif  ! okvan
      
      ! Assign temporary buffers to original arrays
      if (noncolin) then
         ps1_nc(:,:,ibnd) = temp_ps1_nc(:,:)
         ps2_nc(:,:,ibnd,:) = temp_ps2_nc(:,:,:)
      else
         ps1(:,ibnd) = temp_ps1(:)
         ps2(:,ibnd,:) = temp_ps2(:,:)
      end if
   enddo ! nbnd
   !
   !      This term is proportional to beta(k+q+G)
   ! 
   npwq = ngk(ikq)
   if (nkb.gt.0) then
      if (noncolin) then
         call zgemm ('N', 'N', npwq, nbnd*npol, nkb, &
          (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi, npwx)
      else
         call zgemm ('N', 'N', npwq, nbnd, nkb, &
          (1.d0, 0.d0) , vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
      end if
   end if
   !
   !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
   !
   do ikb = 1, nkb
      do ipol = 1, 3
         ok = .false.
         IF (noncolin) THEN
            do ibnd = 1, nbnd
               ok = ok.or.(abs (ps2_nc (ikb, 1, ibnd, ipol) ).gt.eps).or. &
                          (abs (ps2_nc (ikb, 2, ibnd, ipol) ).gt.eps)
            end do
         ELSE
            do ibnd = 1, nbnd
               ok = ok.or. (abs (ps2 (ikb, ibnd, ipol) ) .gt.eps)
            enddo
         ENDIF
         if (ok) then
            do ig = 1, npwq
               igg = igk_k (ig,ikq)
               aux (ig) =  vkb(ig, ikb) * (xk(ipol, ikq) + g(ipol, igg) )
            enddo
            do ibnd = 1, nbnd
               IF (noncolin) THEN
                  call zaxpy(npwq,ps2_nc(ikb,1,ibnd,ipol),aux,1,dvpsi(1,ibnd),1)
                  call zaxpy(npwq,ps2_nc(ikb,2,ibnd,ipol),aux,1, &
                                                          dvpsi(1+npwx,ibnd),1)
               ELSE
                  call zaxpy (npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1)
               END IF
            enddo
         endif
      enddo

   enddo
   deallocate (aux)
   IF (noncolin) THEN
      deallocate (ps2_nc)
      deallocate (ps1_nc)
      deallocate (deff_nc)
      deallocate (temp_ps1_nc)
      deallocate (temp_ps2_nc)
   ELSE
      deallocate (ps2)
      deallocate (ps1)
      deallocate (deff)
      deallocate (temp_ps1)
      deallocate (temp_ps2)
   END IF

   call stop_clock ('dvqpsi_us_on')
   return
end subroutine dvqpsi_us_only
