#!/bin/bash
for m in $(cat multitags.txt)
#do
#    echo $(./split.sh ${m:0:6}_re.txt ${m:6:6}_re.txt everything) 
#done

do
    grep -f ${m:0:6}_re.txt $1 | grep -f ${m:6:6}_re.txt | cut -d ',' -f 2 > ${m}_lintag1.txt

done
