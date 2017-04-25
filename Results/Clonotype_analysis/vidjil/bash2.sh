#/bin/bash
#for entry in `ls -v clone.fa-*` 
for entry in `ls clone.fa-* | sed 's/^\([^0-9]*\)\([0-9]*\)/\1 \2/' | sort -k2,2n | tr -d ' ' |`
do 
	echo "$entry" 
	python readvidjilcloneout.py -f $entry 
done

