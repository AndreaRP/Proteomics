#!/bin/bash
set -e 

# $1 = Directory where the files are
# $2 = File with DE list
directory="${1}/*_formated.txt"
for list in $directory
do
	list_file=$(echo $list | sed 's/_formated.txt/_de_list.txt/')
	> $list_file
	comm -12 <(cut -f1 $list | sort -k1) <(cut $2 -f1 | sort -k1)  >> $list_file
	echo "$list_file done"
done