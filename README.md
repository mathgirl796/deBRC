首先想办法在自己的环境中编译一个KMC3的静态库，因为官网的有时候用不了
```
cd path/to/deBRC
mkdir lib
mv path/to/libkmc_core.a lib/
make
```

使用脚本范例（walk和restore可以添加--useKmerFormat以使用kmerFormat，即每个brc处会输出一整个kmer）
```
#!/bin/bash
#SBATCH -J hg38
#SBATCH -o /home/user/duanran/repo/deBRC/deBRC/experiment/single_genomes/hg38/log/%j.log
#SBATCH -p fat
#SBATCH --mincpus 24
#SBATCH --mem 64G

K=21
USEKMERFORMAT=1
ISMER=0
BASE_DIR=/home/user/duanran/repo/deBRC/deBRC/experiment/single_genomes/hg38
BASE_OUTPUT_NAME=hg38
BASE_INPUT_NAME=hg38.fa
MINCPUS=24
MEM=64

OUTPUT_DIR=$BASE_DIR/k$K
TMP_DIR=$BASE_DIR/k$K/tmp
mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR

/home/user/duanran/repo/deBRC/deBRC/bin/main kmc    -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -k $(expr $K + 1 $(if [ "$USMER" == "1" ];then echo "+ 1";fi)) -w $TMP_DIR $BASE_DIR/$BASE_INPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main convert                    -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main sort   -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
/home/user/duanran/repo/deBRC/deBRC/bin/main split  -t $MINCPUS -r $MEM $(if [ "$ISMER" == "1" ];then echo "--ismer";fi) -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer
if [ "$ISMER" != "1" ]
then
    /home/user/duanran/repo/deBRC/deBRC/bin/main walk   -t $MINCPUS -r $MEM $(if [ "$USEKMERFORMAT" == "1" ];then echo --useKmerFormat;fi) -s asdf -l $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
    # /home/user/duanran/repo/deBRC/deBRC/bin/main restore -t $MINCPUS $(if [ "$USEKMERFORMAT" == "1" ];then echo --useKmerFormat;fi) -s $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME$(if [ "$USEKMERFORMAT" == "1" ];then echo ".kmer";fi).brc
else
    /home/user/duanran/repo/deBRC/deBRC/bin/main unitig -t $MINCPUS -r $MEM -i $OUTPUT_DIR/$BASE_OUTPUT_NAME.i.smer -k $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
fi
rm -f $OUTPUT_DIR/$BASE_OUTPUT_NAME.kmc_pre $OUTPUT_DIR/$BASE_OUTPUT_NAME.kmc_suf
rm -f $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
rm -rf $TMP_DIR
```

注意，restore需要把输入输出和kmer都缓存在内存中，可能需要超多内存（restore实测HG002用了100G，HG002-kmerFormat用了110G，实际可能更多，建议多加几个G）