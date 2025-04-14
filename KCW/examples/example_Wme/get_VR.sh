# Exctract info from the output files of the bare and screened coupomb calculation 
# in a format suitable for easy plotting: V(R)

is=`grep mp1 Si.kcw-screen.in | awk '{print -$3/2}'`
ie=`grep mp1 Si.kcw-screen.in | awk '{print -$3/2+$3-1}'`
js=`grep mp2 Si.kcw-screen.in | awk '{print -$3/2}'`
je=`grep mp2 Si.kcw-screen.in | awk '{print -$3/2+$3-1}'`
ks=`grep mp3 Si.kcw-screen.in | awk '{print -$3/2}'`
ke=`grep mp3 Si.kcw-screen.in | awk '{print -$3/2+$3-1}'`
echo $is $ie
echo $js $je
echo $ks $ke

rm -fr R.dat 
for i in `seq $is $ie`; do 
 for j in `seq $js $je`; do 
  for k in `seq $ks $ke`; do 
   echo $i $j $k  >> R.dat
  done
 done
done

# <ij|V|kl> = \int dr1 dr2 i(r1)*j(r2) v(r1,r2) k(r2)l(r1)
# at the moment only <ij|V|ji> = <rho_i|V|rho_j> = int dr1 dr2 \rho(r1) v(r1,r2) \rho(r2) 
for j in `seq 1 8`; do 
grep "  1           $j           1" barecoulomb.txt  |awk '{print $4, $5}' | sed -r '/^\s*$/d' > pp
paste R.dat pp | awk '{printf "%12.8f %21.18f %21.18f %5i %5i %5i\n", sqrt($1*$1+$2*$2+$3*$3), $4, $5, $1, $2, $3}'> V1${j}bare.dat
sort -n V1${j}bare.dat > pp 
mv pp V1${j}bare.dat
sed -i '1 i\#   |R|             Re[V]                Im[V]              R1    R2    R3' V1${j}bare.dat

grep "  1           $j           1" screencoulomb.txt  |awk '{print $4, $5}' | sed -r '/^\s*$/d' > pp
paste R.dat pp | awk '{printf "%12.8f %21.18f %21.18f %5i %5i %5i\n", sqrt($1*$1+$2*$2+$3*$3), $4, $5, $1, $2, $3}'> V1${j}screen.dat
sort -n V1${j}screen.dat > pp 
mv pp V1${j}screen.dat
sed -i '1 i\#   |R|             Re[W]                Im[W]              R1    R2    R3' V1${j}screen.dat

done 
rm R.dat 

