file=${1:? "usage: file"}
< "$file" grep RI -A 3 | egrep '0|1' > ri.avg.txt
< "$file" grep AR -A 3 | egrep '0|1' > ar.avg.txt
< "$file" grep MI -A 3 | egrep '0|1' > mi.avg.txt
< "$file" grep HI -A 3 | egrep '0|1' > hi.avg.txt

wc -l ri.avg.txt
wc -l ar.avg.txt
wc -l mi.avg.txt
wc -l hi.avg.txt
