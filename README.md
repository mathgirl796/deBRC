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
#SBATCH -J chm13
#SBATCH -o /home/user/duanran/repo/deBRC/deBRC/experiment/chm13/log/%j.out
#SBATCH -p q01
#SBATCH --mincpus=8
#SBATCH --mem=32G
K=31 # 最大为31
BASE_DIR=/home/user/duanran/repo/deBRC/deBRC/experiment/chm13
BASE_OUTPUT_NAME=chm13
BASE_INPUT_NAME=chm13.fa
OUTPUT_DIR=$BASE_DIR/k$K
TMP_DIR=$BASE_DIR/k$K/tmp
mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR
/home/user/duanran/repo/deBRC/deBRC/bin/main kmc -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -k $(expr $K + 1) -r 16 -t 8 -w $TMP_DIR $BASE_DIR/$BASE_INPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main convert -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main sort -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -r 32 -t 8 -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.mer
/home/user/duanran/repo/deBRC/deBRC/bin/main split -o $OUTPUT_DIR/$BASE_OUTPUT_NAME -r 32 -t 8 -w $TMP_DIR $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer
/home/user/duanran/repo/deBRC/deBRC/bin/main walk -s asdf -l $OUTPUT_DIR/$BASE_OUTPUT_NAME.o.smer -t 8 -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $BASE_DIR/$BASE_INPUT_NAME
/home/user/duanran/repo/deBRC/deBRC/bin/main restore -s $OUTPUT_DIR/$BASE_OUTPUT_NAME.smer -t 8 -o $OUTPUT_DIR/$BASE_OUTPUT_NAME $OUTPUT_DIR/$BASE_OUTPUT_NAME.brc
rm -rf $TMP_DIR
```