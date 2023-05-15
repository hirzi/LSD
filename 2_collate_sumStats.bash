#!/bin/bash
dir=$1

find $dir -name "simulatedSumStats*" | sort > $dir/sumStats_filenames.txt

first=$(head -n1 $dir/sumStats_filenames.txt)
head -n1 $first > $dir/sumStats_collated.txt
while read f;do
cat $f | tail -n+2 >> $dir/$1.sumStats_collated.txt
done < $dir/sumStats_filenames.txt