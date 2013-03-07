#!/bin/sh
#
# $Id$

flist="CX_CH_AUX_MAARBLE_ULF_PC12.ceh \
CX_CH_AUX_MAARBLE_ULF_PC35.ceh"


for f in $flist
do
   CLI=1
   while test $CLI -le 4
   do
      fout=`echo $f| sed -e "s=X=${CLI}="`
      echo "Writing $fout"
      cat $f|sed -e "s=XXX=${CLI}=" > $fout
      CLI=$(($CLI+1))
   done
done