# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 23:58:46 2024

@author: zheng

Chapter02

第3个脚本

利用 BioPython 对单个fasta序列读写与常规操作
参考：  Bioinformatics_with_Python_Cookbook Page: P47-49 （并不完全一致）
"""

from Bio import Entrez, SeqIO, Seq

Entrez.email = "zhengshimao007@163.com"

hdl = Entrez.efetch(db='nucleotide', id=['NM_002299'], rettype='fasta')  # Lactase gene

seq = SeqIO.read(hdl, 'fasta')

print(seq)
# ID: NM_002299.4
# Name: NM_002299.4
# Description: NM_002299.4 Homo sapiens lactase (LCT), mRNA
# Number of features: 0
# Seq('AACAGTTCCTAGAAAATGGAGCTGTCTTGGCATGTAGTCTTTATTGCCCTGCTA...GTC')

print(type(seq))
# <class 'Bio.SeqRecord.SeqRecord'>


print(seq.name)
# NM_002299.4
print(seq.description)
# NM_002299.4 Homo sapiens lactase (LCT), mRNA
print(seq.id)
# NM_002299.4
print(seq.seq)
# Seq('AACAGTTCCTAGAAAATGGAGCTGTCTTGGCATGTAGTCTTTATTGCCCTGCTA...GTC')
print(seq.seq.lower)

# 写出单个fasta格式序列
w_hdl = open('example.fasta', 'w')
SeqIO.write(seq, w_hdl, 'fasta')
w_hdl.close()

# 注意：
# 这里SeqIO.write 接受要写入的序列列表（而不仅仅是单个序列）。要小心这个习惯。
# 如果你想写很多序列，不要使用列表list（如前面的代码），因为这将分配大量内存。
# 要么使用迭代器iterator，要么多次使用SeqIO.write函数处理写入的序列。


# 读入单个fasta格式序列
recs = SeqIO.parse('example.fasta', 'fasta')
for rec in recs:
    seq = rec.seq
    print(rec.description)
    print(seq[:10])
    print(seq.alphabet) # 提示没有这个属性了。应该是BioPython版本原因


# 其它方法
part_seq = seq[:10]
print(part_seq)
# AACAGTTCCT

## DNA序列转RNA
rna = part_seq.transcribe()
print(rna)
# AACAGUUCCU

## DNA转小写
lower_seq = part_seq.lower()
print(lower_seq)
# aacagttcct

## DNA小写转大写
up_seq = lower_seq.upper()
print(up_seq)

## 互补序列
comp_seq = part_seq.complement()
print(comp_seq)
# TTGTCAAGGA

## 反向互补序列
rev_comp_seq = part_seq.reverse_complement()
print(rev_comp_seq)
# AGGAACTGTT

## 反向互补序列的RNA
rev_comp_rna = part_seq.reverse_complement_rna()
print(rev_comp_rna)
# AGGAACUGUU

## 翻译成氨基酸 
### 使用默认密码子表
aa = part_seq.translate()
print(aa)
# NSS
# 原序列 AACAGTTCCT
# 翻译分割 AAC AGT TCC T
# 密码子对应的氨基酸 N S S
# (Asn/N)天冬酰胺
# (Ser/S) 丝氨酸

## 统计碱基C的含量
num_c = part_seq.count(sub = 'C')
print(num_c)
# 3

num_g = part_seq.count(sub = 'G')
print(num_g)
# 1

## 统计AC框的数量
num_ac = part_seq.count(sub = 'AC')
print(num_ac) 
# 1

from Bio.Seq import Seq
test = Seq('ACACACA')
print(type(test))
# <class 'Bio.Seq.Seq'>

# 体会两种框的移动方式
test_aca = test.count(sub = 'ACA')
print(test_aca)
## 结果为2，注意体会ACA框的移动

test_aca1 = test.count_overlap(sub = 'ACA')
print(test_aca1)
## 结果为3， 注意体会ACA框的移动

part_seq.endswith('CCT') # 判断序列结尾，返回布尔值。

# True
