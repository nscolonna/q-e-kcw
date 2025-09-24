SUBROUTINE kcw_deallocate_symmetry_arrays()
  !
  ! deallocate arrays related to calculation of symmetries
  !
  USE control_kcw,     ONLY: r, s_w, ft_w, &
                             nsym_w_k, nsym_w_q, nqstot_ibz
  USE control_kcw,     ONLY: xq_ibz, wq_ibz, ibz2fbz, fbz2ibz
  USE control_kcw,     ONLY: nqstar, qstar_iq, qstar_isym, qstar_Gvec
  !
  IMPLICIT NONE
  !
  if( allocated(r)) deallocate(r)
  if( allocated(s_w)) deallocate(s_w)
  if( allocated(ft_w)) deallocate(ft_w)
  if( allocated(nsym_w_k)) deallocate(nsym_w_k)
  if( allocated(nsym_w_q)) deallocate(nsym_w_q)
  if( allocated(nqstot_ibz)) deallocate(nqstot_ibz)
  if( allocated(xq_ibz)) deallocate(xq_ibz)
  if( allocated(wq_ibz)) deallocate(wq_ibz)
  if( allocated(ibz2fbz)) deallocate(ibz2fbz)
  if( allocated(fbz2ibz)) deallocate(fbz2ibz)
  if( allocated ( nqstar )) deallocate(nqstar)
  if( allocated ( qstar_iq )) deallocate(qstar_iq)
  if( allocated ( qstar_isym )) deallocate(qstar_isym)
  if( allocated ( qstar_Gvec )) deallocate(qstar_Gvec)
    !
END SUBROUTINE kcw_deallocate_symmetry_arrays
  
