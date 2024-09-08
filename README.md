首先想办法在自己的环境中编译一个KMC3的静态库，因为官网的有时候用不了
```
cd path/to/deBRC
mkdir lib
mv path/to/libkmc_core.a lib/
make
```

使用脚本范例（walk和restore可以添加--useKmerFormat以使用kmerFormat，即每个brc处会输出一整个kmer）
```
see script/example.sh
```

注意，restore需要把输入输出和kmer都缓存在内存中，可能需要超多内存（restore实测HG002用了100G，HG002-kmerFormat用了110G，实际可能更多，建议多加几个G）