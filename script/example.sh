#!/bin/bash
#SBATCH -J hg38
#SBATCH -o /home/user/duanran/repo/deBRC/deBRC/experiment/single_genomes/hg38/log/%j.log
#SBATCH -p fat
#SBATCH --mincpus 24
#SBATCH --mem 64G

K=21
BRC=1
BRC_USE_KMERFORMAT=1
BRC_USE_IO=1
UNIPATH=0
BASE_DIR=/home/user/duanran/repo/deBRC/deBRC/experiment/single_genomes/hg38
BASE_OUTPUT_NAME=hg38
BASE_INPUT_NAME=hg38.fa
MINCPUS=24
MEM=64

OUTPUT_DIR=$BASE_DIR/k$K
TMP_DIR=$BASE_DIR/k$K/tmp
mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR

/home/user/duanran/repo/deBRC/deBRC/bin/main kmc    -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -k $(expr $K + 1) -w $TMP_DIR $BASE_DIR/$BASE_INPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main convert                    -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main sort   -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
if [ "$BRC" == "1" ];then
    if [ "$BRC_USE_IO" == "1" ];then 
        /home/user/duanran/repo/deBRC/deBRC/bin/main split  -t $MINCPUS -r $MEM --ismer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer
        /home/user/duanran/repo/deBRC/deBRC/bin/main merge  -o $OUTPUT_DIR/$BASE_OUTPUT_NAME.io.smer $OUTPUT_DIR/$BASE_OUTPUT_NAME.i.smer $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer
        /home/user/duanran/repo/deBRC/deBRC/bin/main walk   -t $MINCPUS -r $MEM $(if [ "$BRC_USE_KMERFORMAT" == "1" ];then echo "--useKmerFormat";fi) -s asdf -l $OUTPUT_DIR/$BASE_OUTPUT_NAME.io.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
    else
        /home/user/duanran/repo/deBRC/deBRC/bin/main split  -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer
        /home/user/duanran/repo/deBRC/deBRC/bin/main walk   -t $MINCPUS -r $MEM $(if [ "$BRC_USE_KMERFORMAT" == "1" ];then echo "--useKmerFormat";fi) -s asdf -l $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
    fi
fi
if [ "$UNIPATH" == "1" ];then
    /home/user/duanran/repo/deBRC/deBRC/bin/main split  -t $MINCPUS -r $MEM --ismer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer
    /home/user/duanran/repo/deBRC/deBRC/bin/main unitig -t $MINCPUS -r $MEM -i $OUTPUT_DIR/$BASE_OUTPUT_NAME.i.smer -k $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
fi

rm -f $OUTPUT_DIR/$BASE_OUTPUT_NAME.kmc_pre $OUTPUT_DIR/$BASE_OUTPUT_NAME.kmc_suf
rm -f $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
rm -rf $TMP_DIR
