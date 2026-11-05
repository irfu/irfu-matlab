#!/bin/bash -e
# bepic
# make an empty .cdf file from a skeleton .skt file
#
# uses: skeletoncdf
#
# bepi data in: /data/bepicolombo/reformatter_mirror/pwi/cdf/EFD/L1_prime/202X/
#
# Example: !skeletontable2cdfempty bc_mmo_pwi-efd_l1p_l-e_20251016_r01-v00-00.skt

if [ $# -eq 0 ]; then
    echo " Usage: `basename $0` (sktfilename.skt)"
    exit
fi

wext=cdf

fnnoext=${1%.*}

echo "$fnnoext -> .$wext"

# check for outfilename
if [ -f "${fnnoext}.$wext" ]; then
    mv -v "${fnnoext}.$wext" "$fnnoext.${wext}prev"
fi

skeletoncdf "$fnnoext"

mv  "$fnnoext".cdf "$fnnoext".cdfempty

echo Created file: "$fnnoext".cdfempty

echo Done!
