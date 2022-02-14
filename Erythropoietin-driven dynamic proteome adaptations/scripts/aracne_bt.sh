#!/bin/bash
tfFile="../data/tf.txt"
inFolder="../aracne_bt/inputs/"
outFolder="../aracne_bt/outputs/"

N=4
cd $inFolder
for f in *.txt;do
        java -jar /home/gandri/Programs/aracne2.jar -i ${inFolder}$f -o ${outFolder}$f -e 0 -t 0.35 -l ${tfFile}&
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge $N ]; then
                wait
                NPROC=0
        fi
done

# Rename .txt 2 .adj
cd $outFolder
for file in *.txt; do
    mv "$file" "`basename $file .txt`.adj"
done