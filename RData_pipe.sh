#!/bin/bash

REF="/home/UTHSCSA/hsiaot/src/RSEM_ref/Ensembl_Homo_sapiens_69"
echo RSEM Reference Path = $REF
INPUTPREFIX="/home/UTHSCSA/hsiaot/PECK/RNA_pathseq/"
OUTPUTPREFIX="RSEM"
# Check OUTPUT directory exists
if [ ! -d "$OUTPUTPREFIX" ]; then
   # Control will enter here if $DIRECTORY doesn't exist
  mkdir -m 777 RSEM/
fi
echo Input Path = $INPUTPREFIX
echo Output Path = $OUTPUTPREFIX

SAMPLES=($(find $INPUTPREFIX -type d -name 'ERR*'))
ELEMENTS=${#SAMPLES[@]}

# for loop
for (( i=0;i<$ELEMENTS;i++)); do
   echo ${SAMPLES[$i]}
   SAMPLE=($(basename ${SAMPLES[$i]}))
   echo $SAMPLE
   OUTPUTDIR="$OUTPUTPREFIX/$SAMPLE"
   mkdir $OUTPUTDIR
   R1FILES=($(find ${SAMPLES[$i]} -type f -iname '*_1.fastq.gz'))
   R2FILES=($(find ${SAMPLES[$i]} -type f -iname '*_2.fastq.gz'))
   NUMPAIRS=${#R2FILES[@]}
   echo $NUMPAIRS
   if [ $NUMPAIRS == 0 ]; then
      # single read
      echo "Single Reads"
      NUMREADS=${#R1FILES[@]}
      for ((j=0;j<$NUMREADS;j++)); do
        echo cd $OUTPUTDIR 

