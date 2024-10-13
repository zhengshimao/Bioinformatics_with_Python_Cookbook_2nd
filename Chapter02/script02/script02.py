# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 22:43:23 2024

@author: zheng

Chapter02

第2个脚本

利用 BioPython 下载GenBank格式中的参考文献
参考：  Bioinformatics_with_Python_Cookbook Page: P46 （并不完全一致）

"""

## 注意： 从下方代码看，并不是所有gb都能有pubmed_id记录

from Bio import Entrez, SeqIO

Entrez.email = "zhengshimao007@163.com"

handle = Entrez.esearch(db="nucleotide", term='KM288867')

rec_list = Entrez.read(handle)

hdl = Entrez.efetch(db = 'nucleotide', id = rec_list['IdList'], rettype = "gb")

print(hdl) # <_io.TextIOWrapper encoding='UTF-8'>

rec = list(SeqIO.parse(hdl, format = "gb"))

print(rec)
# [SeqRecord(seq=Seq('ATATGTAAAACCAAAATAAATTAAACAGAATTTATTTTTAAAAGATTTATTTGT...CAT'), id='KM288867.1', name='KM288867', description='Plasmodium falciparum clone PF3D7_0709000 chloroquine resistance transporter (CRT) gene, complete cds', dbxrefs=[])]
print(type(rec)) # <class 'list'>

print(type(rec[0])) # <class 'Bio.SeqRecord.SeqRecord'>
print(rec[0])
# ID: KM288867.1
# Name: KM288867
# Description: Plasmodium falciparum clone PF3D7_0709000 chloroquine resistance transporter (CRT) gene, complete cds
# Number of features: 23
# /molecule_type=DNA
# /topology=linear
# /data_file_division=INV
# /date=12-NOV-2014
# /accessions=['KM288867']
# /sequence_version=1
# /keywords=['']
# /source=Plasmodium falciparum (malaria parasite P. falciparum)
# /organism=Plasmodium falciparum
# /taxonomy=['Eukaryota', 'Sar', 'Alveolata', 'Apicomplexa', 'Aconoidasida', 'Haemosporida', 'Plasmodiidae', 'Plasmodium', 'Plasmodium (Laverania)']
# /references=[Reference(title='Versatile control of Plasmodium falciparum gene expression with an inducible protein-RNA interaction', ...), Reference(title='Direct Submission', ...)]
# Seq('ATATGTAAAACCAAAATAAATTAAACAGAATTTATTTTTAAAAGATTTATTTGT...CAT')


from Bio import Medline

refs = rec[0].annotations['references']

f_hdl = open("KM288867_PubMed_infor.txt","w+")

for ref in refs:
    if ref.pubmed_id != '': # 从这个代码看，并不是所有gb都能有pubmed_id记录
        print(ref.pubmed_id)
        handle = Entrez.efetch(db = 'pubmed', id = [ref.pubmed_id], 
                               rettype = 'medline', retmode = 'text')
        records = Medline.parse(handle)
        # print(type(records)) # generator 生成器 不是很理解
        # generator保存的是算法，每次调用next()，就计算出下一个元素的值，
        # 直到计算到最后一个元素，没有更多的元素时，抛出StopIteration的错误。
        for med_rec in records:
            for k,v in med_rec.items():
                print('%s: %s' % (k,v))
                print('%s: %s' % (k,v), file= f_hdl)
        f_hdl.close()
# PMID: 25370483
# OWN: NLM
# STAT: MEDLINE
# DCOM: 20160112
# LR: 20181113
# IS: 2041-1723 (Electronic) 2041-1723 (Linking)
## VI 是文章的卷 volume 
# VI: 5
## 如果有期（issue），则用IP 表示
## 文章发表日期
# DP: 2014 Nov 5
## TI为文章DOI
# TI: Versatile control of Plasmodium falciparum gene expression with an inducible protein-RNA interaction.
## PG 是文章的页码，对于NC是指文章编号（Article number）
# PG: 5329
## LID为文章DOI
# LID: 10.1038/ncomms6329 [doi] 
# AB: The available tools for conditional gene expression in Plasmodium falciparum are limited. Here, to enable reliable control of target gene expression, we build a system to efficiently modulate translation. We overcame several problems associated with other approaches for regulating gene expression in P. falciparum. Specifically, our system functions predictably across several native and engineered promoter contexts, and affords control over reporter and native parasite proteins irrespective of their subcellular compartmentalization. Induction and repression of gene expression are rapid, homogeneous and stable over prolonged periods. To demonstrate practical application of our system, we used it to reveal direct links between antimalarial drugs and their native parasite molecular target. This is an important outcome given the rapid spread of resistance, and intensified efforts to efficiently discover and optimize new antimalarial drugs. Overall, the studies presented highlight the utility of our system for broadly controlling gene expression and performing functional genetics in P. falciparum.
## FAU 为文章作者较全的写法
# FAU: ['Goldfless, Stephen J', 'Wagner, Jeffrey C', 'Niles, Jacquin C']
## AU 为文章作者简写
# AU: ['Goldfless SJ', 'Wagner JC', 'Niles JC']
## AD 为文章作者所在单位
# AD: ['Department of Biological Engineering, Massachusetts Institute of Technology, 77 Massachusetts Avenue, Cambridge, Massachusetts 02139, USA.', 'Department of Biological Engineering, Massachusetts Institute of Technology, 77 Massachusetts Avenue, Cambridge, Massachusetts 02139, USA.', 'Department of Biological Engineering, Massachusetts Institute of Technology, 77 Massachusetts Avenue, Cambridge, Massachusetts 02139, USA.']
# LA: ['eng']
# SI: ['GENBANK/KM288848', 'GENBANK/KM288849', 'GENBANK/KM288850', 'GENBANK/KM288851', 'GENBANK/KM288852', 'GENBANK/KM288853', 'GENBANK/KM288854', 'GENBANK/KM288855', 'GENBANK/KM288856', 'GENBANK/KM288857', 'GENBANK/KM288858', 'GENBANK/KM288859', 'GENBANK/KM288860', 'GENBANK/KM288861', 'GENBANK/KM288862', 'GENBANK/KM288863', 'GENBANK/KM288864', 'GENBANK/KM288865', 'GENBANK/KM288866', 'GENBANK/KM288867']
# GR: ['1DP2OD007124/OD/NIH HHS/United States', 'DP2 OD007124/OD/NIH HHS/United States', '5-T32-ES007020/ES/NIEHS NIH HHS/United States', 'T32 GM008334/GM/NIGMS NIH HHS/United States', 'T32 ES007020/ES/NIEHS NIH HHS/United States', 'P30 ES002109/ES/NIEHS NIH HHS/United States', '5-T32-GM08334/GM/NIGMS NIH HHS/United States']
## PT应该是指文章类型
# PT: ['Journal Article', 'Research Support, N.I.H., Extramural', "Research Support, Non-U.S. Gov't"]
# DEP: 20141105
# PL: England
## 发表期刊简写
# TA: Nat Commun
## 发表期刊全称
# JT: Nature communications
# JID: 101528555
# RN: ['0 (Aptamers, Nucleotide)']
# SB: IM
# MH: ['Aptamers, Nucleotide', 'Base Sequence', '*Gene Expression Regulation', '*Genetic Techniques', 'Molecular Sequence Data', 'Plasmodium falciparum/genetics/*metabolism']
## 文章在Pubmed Central的编号
# PMC: PMC4223869
# MID: ['NIHMS630149']
# COIS: ['COMPETING FINANCIAL INTERESTS S.J.G and J.C.N. are inventors of the genetically', 'encoded protein-binding RNA aptamer technology described and have filed patent', 'applications with other co-inventors.']
# EDAT: 2014/11/06 06:00
# MHDA: 2016/01/13 06:00
# PMCR: ['2015/05/05']
# CRDT: ['2014/11/06 06:00']
# PHST: ['2014/04/15 00:00 [received]', '2014/09/20 00:00 [accepted]', '2014/11/06 06:00 [entrez]', '2014/11/06 06:00 [pubmed]', '2016/01/13 06:00 [medline]', '2015/05/05 00:00 [pmc-release]']
# AID: ['ncomms6329 [pii]', '10.1038/ncomms6329 [doi]']
# PST: epublish
# SO: Nat Commun. 2014 Nov 5;5:5329. doi: 10.1038/ncomms6329.

## 测试缩写含义
## PMID: 39380220
f_hdl2 = open("PMID_39380220_kiwifruit.txt","w+")
for ref in refs:
    if ref.pubmed_id != '':
        print(ref.pubmed_id)
        handle = Entrez.efetch(db = 'pubmed', id = [39380220], 
                               rettype = 'medline', retmode = 'text')
        records = Medline.parse(handle)
        # print(type(records)) # generator 生成器 不是很理解
        # generator保存的是算法，每次调用next()，就计算出下一个元素的值，
        # 直到计算到最后一个元素，没有更多的元素时，抛出StopIteration的错误。
        for med_rec in records:
            for k,v in med_rec.items():
                print('%s: %s' % (k,v))
                print('%s: %s' % (k,v), file= f_hdl2)
        f_hdl2.close()




