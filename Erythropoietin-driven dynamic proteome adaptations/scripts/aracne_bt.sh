#!/bin/bash
tfFile="../data/tf.txt"
inFolder="../aracne_bt/inputs/"
outFolder="../aracne_bt/outputs/"

N=4
cd $inFolder
for f in *.txt
do
	((i=i%N)); ((i++==0)) && wait
	java -jar /Users/gandrieux/Programs/aracne2.jar -i ${inFolder}$f -o ${outFolder}$f -e 0 -l ${tfFile}&
done	
wait

# Rename .txt 2 .adj
cd $outFolder
for file in *.txt; do
    mv "$file" "`basename $file .txt`.adj"
done