!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE init_nsg
   !-----------------------------------------------------------------------
   !
   ! This routine computes the starting ns (for DFT+U+V calculation) filling
   ! up the Hubbard manifold (we are only interested in the on-site potential 
   ! for the moment) according to the Hund's rule (valid for the isolated atoms 
   ! on which starting potential is built), and to the starting_magnetization:
   ! majority spin levels are populated first, then the remaining electrons
   ! are equally distributed among the minority spin states
   !
   USE kinds,       ONLY : DP
   USE ions_base,   ONLY : nat, ityp
   USE uspp_param,  ONLY : upf
   USE lsda_mod,    ONLY : nspin, starting_magnetization
   USE ldaU,        ONLY : Hubbard_l, Hubbard_l2, Hubbard_l3, hubbard_occ, &
                           nsg, ldim_u, backall, is_hubbard, is_hubbard_back
   !
   IMPLICIT NONE
   REAL(DP) :: totoc, totoc_b
   INTEGER :: ldim, ldim2, na, nt, is, m1, majs, mins, viz
   LOGICAL :: nm        ! true if the atom is non-magnetic
   INTEGER, EXTERNAL :: find_viz
   !
   nsg(:,:,:,:,:) = (0.d0, 0.d0)
   ! 
   DO na = 1, nat
      ! 
      ! The index of atom 'na' in the neighbors list
      !
      viz = find_viz(na,na) 
      !
      nt = ityp(na)
      !
      IF ( is_hubbard(nt) ) THEN 
         !
         ldim = 2*Hubbard_l(nt)+1
         nm = .TRUE.
         !
         totoc = hubbard_occ(nt,1)
         !
         IF (nspin.EQ.2) THEN
            IF (starting_magnetization(nt).GT.0.d0) THEN 
               nm = .FALSE.
               majs = 1  
               mins = 2  
            ELSEIF (starting_magnetization(nt).LT.0.d0) THEN 
               nm = .FALSE.
               majs = 2  
               mins = 1  
            ENDIF  
         ENDIF
         !
         IF (.NOT.nm) THEN
            ! Atom is magnetic 
            IF (totoc.GT.ldim) THEN
               DO m1 = 1, ldim
                  nsg (m1,m1,viz,na,majs) = 1.d0
                  nsg (m1,m1,viz,na,mins) = (totoc - ldim) / ldim
               ENDDO
            ELSE
               DO m1 = 1, ldim
                  nsg (m1,m1,viz,na,majs) = totoc / ldim
               ENDDO
            ENDIF
         ELSE  
            ! Atom is non-magnetic
            DO is = 1, nspin
               DO m1 = 1, ldim  
                  nsg (m1,m1,viz,na,is) = totoc /  2.d0 / ldim
               ENDDO  
            ENDDO  
         ENDIF  
         !
         ! Background part
         !
         IF ( is_hubbard_back(nt) ) THEN
            !     
            IF (.NOT.backall(nt)) THEN
               ! Fill in the second Hubbard manifold
               ldim2 = 2*Hubbard_l2(nt)+1
               totoc_b = hubbard_occ(nt,2)
               DO is = 1, nspin
                  DO m1 = ldim+1, ldim_u(nt)
                     nsg (m1,m1,viz,na,is) = totoc_b / 2.d0 / ldim2
                  ENDDO
               ENDDO
            ELSE
               ! Fill in the second Hubbard manifold
               ldim2 = 2*Hubbard_l2(nt)+1
               totoc_b = hubbard_occ(nt,2)
               DO is = 1, nspin
                  DO m1 = ldim+1, ldim+2*Hubbard_l2(nt)+1
                     nsg (m1,m1,viz,na,is) = totoc_b / 2.d0 / ldim2
                  ENDDO
               ENDDO
               ! Fill in the third Hubbard manifold
               ldim2 = 2*(Hubbard_l2(nt) + Hubbard_l3(nt)) + 2
               totoc_b = hubbard_occ(nt,3)
               DO is = 1, nspin
                  DO m1 = ldim+2*Hubbard_l2(nt)+2, ldim_u(nt)
                     nsg (m1,m1,viz,na,is) = totoc_b / 2.d0 / ldim2
                  ENDDO
               ENDDO
            ENDIF
            !
         ENDIF
         !
      ENDIF
      ! 
   ENDDO ! na 
   !
   RETURN
   !  
END SUBROUTINE init_nsg
!-----------------------------------------------------------------------
