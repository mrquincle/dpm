#!/bin/sh

file=${1:? "usage: file output-path"}
opath=${2:? "usage: file output-path"}

mkdir -p $opath

# old format is 3, new is 1
format=3

< "$file" grep RI -A $format | egrep '0|1' > $opath/ri.avg.txt
< "$file" grep AR -A $format | egrep '0|1' > $opath/ar.avg.txt
< "$file" grep MI -A $format | egrep '0|1' > $opath/mi.avg.txt
< "$file" grep HI -A $format | egrep '0|1' > $opath/hi.avg.txt

wc -l $opath/ri.avg.txt
wc -l $opath/ar.avg.txt
wc -l $opath/mi.avg.txt
wc -l $opath/hi.avg.txt
