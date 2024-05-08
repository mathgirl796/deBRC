首先想办法在自己的环境中编译一个KMC3，因为官网的有时候用不了
```
mkdir lib
mv path/to/libkmc_core.a lib/
make
```


```
#!/bin/bash

# 创建临时文件夹，防止多个任务的临时文件冲突
mkdir -p tmp/test

# 调用KMC3 API，进行kmer counting，k最大为32；这一步会生成文件test.kmc_pre和test.kme_suf（文件名由-o参数加后缀而得）
# 这一步内存控制比较好，因为不是我控制的，是KMC3控制的
# 限制使用16G内存的情况下，在chm13基因组上用了14.286636GB内存，用时1分多
bin/main kmc -o test -k 22 -r 16 -t 8 -w tmp/test test.fa.gz

# 把KMC3格式的kmer转换为.mer文件格式，具体定义为：uint32_t kmerLength | uint64_t kmerNum | uint64_t kmer[kmerNum]
# 这一步会生成文件test.mer（最后的参数test是在指定test.kmc_pre和test.kmc_suf，这是KMC3的API要求）
# 这一步基本不用内存，chm13基因组内存消耗0.055645GB，用时6分多
bin/main convert -o test test

# 这一步把.mer文件进行排序，生成文件test.smer
# 这一步内存似乎有点溢出，限制16G结果用了25.036301GB，用时9分多
bin/main sort -o test -r 16 -t 8 -w tmp/test test.mer

# 这一步查找.smer文件中的多出kmer，生成文件test.o.smer，该文件的k将比输入文件小1。如果指定--kmer，还会输出一个test.k.smer，用于brc的恢复
# 这一步预计的内存用量约等于.o.smer的大小，如果要求输出.k.smer的话还得把它的大小加上
# 因为临时结果直接存在内存里
# chm13只生成.o.smer内存消耗2.790329GB，用时3分多
bin/main split -o test test.smer

# 这一步根据test.fa和test.o.smer，生成test.brc。目前仅支持fasta文件在自己所创建的索引上生成brc，不支持一个fasta在另一个fasta创建的索引上生成brc
# 这一步的-t参数（线程数）建议设置为23，同时处理人的23个主要染色体
# 这步内存消耗约等于fa解压后的大小 + 索引大小 + 结果大小，因为全部加载/缓存进内存了
# chm13上内存消耗7.364403GB，8个worker用时16分多，23个worker的话十来分钟应该就干完了（最长的chr1所用时间）
bin/main walk -k asdf -l test.o.smer -o test -t 8 test.fa.gz

# 删除临时文件夹
rm -rf /home/user/duanran/repo/deBRC/deBRC/experiment/tmp/chm13.draft_v1.1
```

在chm13上的输出如下
```
-rw-rw-r--. 1 duanran duanran 107M May  7 16:00 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.brc
-rw-rw-r--. 1 duanran duanran 1.1M May  7 12:59 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.kmc_pre
-rw-rw-r--. 1 duanran duanran  14G May  7 12:59 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.kmc_suf
-rw-rw-r--. 1 duanran duanran  18G May  7 13:04 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.mer
-rw-rw-r--. 1 duanran duanran 473M May  7 13:15 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.o.smer
-rw-rw-r--. 1 duanran duanran 107M May  7 15:44 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.passN.brc
-rw-rw-r--. 1 duanran duanran  18G May  7 13:12 /home/user/duanran/repo/deBRC/deBRC/experiment/output/chm13.draft_v1.1.smer
```
其中，.brc和.passN.brc经过diff命令检查，完全一致