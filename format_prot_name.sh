#!/bin/bash
set -e 

# $1 = Directory where the files are
directory="${1}/*.csv"
for list in $directory
do
	formated_file=$(echo $list | sed 's/.csv/_formated.txt/')
	> $formated_file
	cat $list | sed '1 s/Row Labels/ID1\tID2/' | sed 's/_[^\t]*//' | sed 's/>sp|//g' | sed 's/|/\t/g' >> $formated_file
	#| cut -d $'\t' -f2,4,5,6,7,8,9 
	echo "$formated_file done"
done
