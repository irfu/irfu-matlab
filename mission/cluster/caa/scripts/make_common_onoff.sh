#!/bin/sh

# Make a common file from COVERAGE_##_01.TXT
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$

sc_list="1 \
2 \
3 \
4"
tmp_f=common_tmp.txt

rm -rf $tmp_f; touch $tmp_f

for cli in $sc_list
do
	cat EFWONOFF_$cli'.TXT' | grep EFWon > on_tmp
	cat EFWONOFF_$cli'.TXT' | grep EFWoff > off_tmp
	cat on_tmp |awk "{print \$1 \"T\" \$2 \".000Z $cli 1\"}" >> $tmp_f
	cat off_tmp|awk "{print \$1 \"T\" \$2 \".000Z $cli 0\"}" >> $tmp_f
	rm -rf on_tmp off_tmp
done

cat $tmp_f| sort > EFWONOFF_COMM.dat
rm $tmp_f
