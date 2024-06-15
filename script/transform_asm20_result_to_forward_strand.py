# 说明：首先用 minimap2 -ax asm5 把样本基因组比对到chm13上得到sam文件，再用该脚本从该sam文件还原出全部比对到正链上的样本基因组

import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm
import sys

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <alignment_file> <output_file>")

bf = pysam.AlignmentFile(sys.argv[1], 'rb')

total_seq_num = 0
total_base_num = 0
rc_seq_num = 0
records = []
for r in tqdm(bf):
    if not r.is_secondary and not r.is_supplementary:
        total_seq_num += 1
        total_base_num += len(r.query_sequence)
        if r.is_reverse:
            info = f"{r.reference_name}, {r.query_name}"
            # print(info, end="\b" * len(info))
            rc_seq_num += 1
            records.append(SeqRecord(Seq(r.query_sequence), id=r.query_name, description="")) # SAMv1.pdf:"Bit 0x10 indicates whether SEQ has been reverse complemented and QUAL reversed. "
        else:
            records.append(SeqRecord(Seq(r.query_sequence), id=r.query_name, description=""))

print(f"total_seq_num: {total_seq_num}, total_base_num: {total_base_num}, rc_seq_num: {rc_seq_num}")
SeqIO.write(records, sys.argv[2], "fasta")

bf.close()