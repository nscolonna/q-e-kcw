#!/bin/bash
# compare results/ with reference/ for example_2D (h-BN monolayer, 2D cutoff)
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

tol_ene=0.000001
tol_eig=0.001
tol_eps=0.0001
tol_r2=0.00001
tol_a=0.0001
echo
echo " Testing the results ..."
echo "  tolerances: "
printf "     energy %10.7f Ry  \n" $tol_ene
printf "     eigval %10.7f eV  \n" $tol_eig
printf "     eps    %10.7f     \n" $tol_eps
printf "     R^2    %10.7f A^2 \n" $tol_r2
printf "     alpha  %10.7f     \n" $tol_a

# cmp label a a_ref tol : prints a warning and sets check=0 if |a-a_ref| > tol
cmp() {
  err=`echo $2 $3 | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
  if (( $(echo "$err > $4" | bc -l) )); then
    echo -e "${RED}  WARNING: ${NC} $1 does not match: $2 $3 $err"
    check=0
  fi
}

## PWSCF (2D cutoff): total energy, VBM, CBM
check=1
cmp "PWSCF total energy" `grep ! results/hBN.scf.out | awk '{print $5}'` `grep ! reference/hBN.scf.out | awk '{print $5}'` $tol_ene
cmp "PWSCF VBM" `grep high results/hBN.scf.out | awk '{print $7}'` `grep high reference/hBN.scf.out | awk '{print $7}'` $tol_eig
cmp "PWSCF CBM" `grep high results/hBN.scf.out | awk '{print $8}'` `grep high reference/hBN.scf.out | awk '{print $8}'` $tol_eig
if (( $check )); then echo -e "  ${GREEN}PWSCF       OK!${NC}"; fi

## PHONON: dielectric tensor (diagonal)
check=1
for i in 1 2 3; do
  cmp "eps($i,$i)" `sed -n "${i}p" results/eps.dat | awk -v i=$i '{print $i}'` `sed -n "${i}p" reference/eps.dat | awk -v i=$i '{print $i}'` $tol_eps
done
if (( $check )); then echo -e "  ${GREEN}PHONON      OK!${NC}"; fi

## W90: spreads of the 4 occupied and 2 empty Wannier functions
check=1
for i in `seq 2 5`; do
  cmp "R^2 occ orbital $((i-1))" `grep -A 4 "Final State" results/wann_occ/hBN.wout | sed -n "${i}p" | awk '{print $11}'` \
                                 `grep -A 4 "Final State" reference/wann_occ/hBN.wout | sed -n "${i}p" | awk '{print $11}'` $tol_r2
done
for i in 2 3; do
  cmp "R^2 emp orbital $((i+3))" `grep -A 2 "Final State" results/wann_emp/hBN_emp.wout | sed -n "${i}p" | awk '{print $11}'` \
                                 `grep -A 2 "Final State" reference/wann_emp/hBN_emp.wout | sed -n "${i}p" | awk '{print $11}'` $tol_r2
done
if (( $check )); then echo -e "  ${GREEN}WANNIER     OK!${NC}"; fi

## KC_SCREEN: q+G=0 terms of the 2D scheme and alpha of the 6 orbitals
check=1
for s in "uPi(q+G=0) estimation" "rPi(q+G=0) estimation"; do
  cmp "$s" `grep -F "$s" results/hBN.kcw-screen.out | awk '{print $4}'` `grep -F "$s" reference/hBN.kcw-screen.out | awk '{print $4}'` $tol_a
done
for i in `seq 1 6`; do
  cmp "alpha orbital $i" `grep "relaxed =" results/hBN.kcw-screen.out | sed -n "${i}p" | awk '{print $12}'` \
                         `grep "relaxed =" reference/hBN.kcw-screen.out | sed -n "${i}p" | awk '{print $12}'` $tol_a
done
if (( $check )); then echo -e "  ${GREEN}KC_SCREEN   OK!${NC}"; fi

## KC_HAM: KI eigenvalues at the k points of the 3x3x1 mesh
check=1
nk=`grep -c "KI  " reference/hBN.kcw-ham.out`
for i in `seq 1 $nk`; do
  a=(`grep "KI  " results/hBN.kcw-ham.out | sed -n "${i}p"`)
  b=(`grep "KI  " reference/hBN.kcw-ham.out | sed -n "${i}p"`)
  for j in `seq 1 $(( ${#b[@]} - 1 ))`; do
    cmp "KI eig ik=$i,v=$j" ${a[$j]} ${b[$j]} $tol_eig
  done
done
if (( $check )); then echo -e "  ${GREEN}KC_HAM      OK!${NC}"; fi
echo
