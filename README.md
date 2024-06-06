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

if [ $# -lt 7 ]; then
    echo "Usage: $(basename $0) <K> <USE_KMER_FORMAT(1 for yes)> <BASE_DIR> <BASE_OUTPUT_NAME(will output to \$BASE_NAME/k\$K/\$BASE_OUPTUT_NAME.*)> <BASE_INPUT_NAME(will reference to \$BASE_NAME/\$BASE_INPUT_NAME)> <MINCPUS> <MEM(GB)>"
    exit 1
fi
echo "shell parameter (num: $#): [$1][$2][$3][$4][$5][$6][$7]"
K=$1
USEKMERFORMAT=$2
BASE_DIR=$3
BASE_OUTPUT_NAME=$4
BASE_INPUT_NAME=$5
MINCPUS=$6
MEM=$7

OUTPUT_DIR=$BASE_DIR/k$K
TMP_DIR=$BASE_DIR/k$K/tmp
mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR
/home/user/duanran/repo/deBRC/deBRC/bin/main kmc    -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -k $(expr $K + 1) -w $TMP_DIR $BASE_DIR/$BASE_INPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main convert           -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME
rm -f $OUTPUT_DIR/$BASE_OUTPUT_NAME.kmc_pre $OUTPUT_DIR/$BASE_OUTPUT_NAME.kmc_suf
/home/user/duanran/repo/deBRC/deBRC/bin/main sort   -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
rm -f $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
/home/user/duanran/repo/deBRC/deBRC/bin/main split  -t $MINCPUS -r $MEM -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer
/home/user/duanran/repo/deBRC/deBRC/bin/main walk   -t $MINCPUS -r $MEM $(if [ "$2" -eq 1 ];then echo --useKmerFormat;fi) -s asdf -l $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
# /home/user/duanran/repo/deBRC/deBRC/bin/main restore -t $MINCPUS $(if [ "$2" -eq 1 ];then echo --useKmerFormat;fi) -s $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME$(if [ "$2" -eq 1 ];then echo .kmer;fi).brc
rm -rf $TMP_DIR
```
注意，restore需要把输入输出和kmer都缓存在内存中，可能需要超多内存（restore实测HG002用了100G，HG002-kmerFormat用了110G，实际可能更多，建议多加几个G）