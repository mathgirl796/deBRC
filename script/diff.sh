#!/bin/bash

tr -d '\n' < $1 > $1.diff.tmp
tr -d '\n' < $2 > $2.diff.tmp
echo $(tr -d '\n' < $1) > $1.diff.comb
echo $(tr -d '\n' < $2) >> $1.diff.comb
diff $1.diff.tmp $2.diff.tmp
# rm $1.diff.tmp $2.diff.tmp