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
bin/main kmc -o test -k 22 -r 16 -t 8 -w tmp/test test.fa.gz

# 把KMC3格式的kmer转换为.mer文件格式，具体定义为：uint32_t kmerLength | uint64_t kmerNum | uint64_t kmer * n
# 这一步会生成文件test.mer（最后的参数test是在指定test.kmc_pre和test.kmc_suf，这是KMC3的API要求）
bin/main convert -o test test

# 这一步把.mer文件进行排序，生成文件test.smer
bin/main sort -o test -r 16 -t 8 test.mer

# 这一步查找.smer文件中的多出kmer，生成文件test.o.smer，该文件的k将比输入文件小1。如果指定--kmer，还会输出一个test.k.smer，用于brc的恢复
bin/main split -o test test.smer

# 这一步根据test.fa和test.o.smer，生成test.brc。目前仅支持fasta文件在自己所创建的索引上生成brc，不支持一个fasta在另一个fasta创建的索引上生成brc
# 这一步的-t参数（线程数）建议设置为23，同时处理人的23个主要染色体
bin/main walk -k asdf -l test.o.smer -o test -t 8 test.fa.gz

# 删除临时文件夹
rm -rf /home/user/duanran/repo/deBRC/deBRC/experiment/tmp/chm13.draft_v1.1
```