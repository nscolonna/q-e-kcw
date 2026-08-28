#!/bin/sh
#
# extract_proj.sh -- build a *.proj file (as read by epsilon.x, calculation='jdos',
# via the filproj1/filproj2 namelist variables) out of the standard ASCII
# projection file written by projwfc.x (prefix.projwfc_up / prefix.projwfc_down,
# see PP/src/write_proj.f90::write_proj_file).
#
# For every Kohn-Sham state (ik,ibnd) the output file contains the SUM of
# |<atomic wfc|psi>|^2 over the atomic wavefunctions selected by <atom_list>,
# <orb> and <ml_list>. This is the per-(ik,ibnd) weight used by epsilon.x to
# build the projected joint density of states (pjdos).
#
# <ml_list> selects which magnetic quantum number(s) m of <orb> to sum over
# (the last integer column of the header line, "m" in write_proj_file, 1..2l+1
# in the real spherical harmonics convention used by projwfc.x). Use "all" to
# sum over every m (e.g. the full d shell); use a comma separated subset to
# isolate individual orbitals or symmetry-adapted combinations, e.g. for a
# d shell in Oh symmetry with the projwfc.x convention
#   1 dz2, 2 dzx, 3 dzy, 4 dx2-y2, 5 dxy
# so ml_list="2,3,5" gives t2g and ml_list="1,4" gives eg: this is what makes
# it possible to inspect a specific d-d (e.g. t2g -> eg) transition instead
# of the full, m-summed, d -> d jdos.
#
# USAGE:
#   extract_proj.sh <projwfc_ascii_file> <nk> <nbnd> <atom_list> <orb> <ml_list> <outfile>
#
#     <projwfc_ascii_file>  e.g. LiCoO2.projwfc_up
#     <nk>                  number of k points (same as in the pw.x nscf run)
#     <nbnd>                number of bands (same as in the pw.x nscf run)
#     <atom_list>           comma separated list of atom indices (the "na"
#                            column of the header lines), e.g. "1" or "2,3"
#     <orb>                 orbital label as printed by projwfc.x, e.g. 3D, 2P
#     <ml_list>             comma separated list of m values to sum over,
#                            e.g. "1,4" for eg, "2,3,5" for t2g, or "all"
#     <outfile>             name of the file to write (ik ibnd sum_weight)
#
if [ $# -ne 7 ]; then
    echo "USAGE: extract_proj.sh <projwfc_file> <nk> <nbnd> <atom_list> <orb> <ml_list> <outfile>"
    exit 1
fi

FILE=$1
NK=$2
NBND=$3
ATOMS=$4
ORB=$5
MLIST=$6
OUT=$7

awk -v atoms="$ATOMS" -v orb="$ORB" -v mlist="$MLIST" '
BEGIN {
    n = split(atoms, a, ",")
    for (i = 1; i <= n; i++) want_atom[a[i]] = 1
    all_m = (mlist == "all")
    if (!all_m) {
        nm = split(mlist, ml, ",")
        for (i = 1; i <= nm; i++) want_m[ml[i]+0] = 1
    }
    match_block = 0
    norder = 0
}
{
    # header line of an atomic wfc block:  nwfc  na  atm  els   n   l   m
    if (NF == 7 && $3 ~ /^[A-Za-z]+$/ && $4 ~ /^[0-9][A-Za-z]$/) {
        match_block = (($2 in want_atom) && $4 == orb && (all_m || ($7 in want_m)))
        next
    }
    # data line of an atomic wfc block:  ik  ibnd  |<phi|psi>|^2
    if (match_block && NF == 3 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/) {
        key = $1 "," $2
        sum[key] += $3
        if (!(key in seen)) { seen[key] = 1; order[++norder] = key }
        next
    }
}
END {
    if (norder == 0) {
        print "extract_proj.sh: no matching atomic wfc block found (atoms=" atoms ", orb=" orb ", ml=" mlist ")" > "/dev/stderr"
        exit 1
    }
    for (i = 1; i <= norder; i++) {
        split(order[i], ik_ib, ",")
        printf "%8d%8d%20.10f\n", ik_ib[1], ik_ib[2], sum[order[i]]
    }
}' "$FILE" > "$OUT"

NLINES=$(wc -l < "$OUT")
NEXP=$((NK*NBND))
if [ "$NLINES" -ne "$NEXP" ]; then
    echo "extract_proj.sh: WARNING $OUT has $NLINES lines, expected $NEXP (nk=$NK x nbnd=$NBND)"
fi
