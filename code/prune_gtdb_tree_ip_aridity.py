#!/usr/bin/env python3

from ete4 import Tree
import gzip

#CLARA TREE
with gzip.open("/support/databases/gtdb/genomic_files_reps_v226/aux/bac120.tree.gz", "rt",encoding='utf-8') as f:
    t = Tree(f.read(), parser=1) #GTDB format has internal node labels, needs format=1 to be read
# load IDs to keep
with open("scripts/JAE_HK_Aridity/data/clara_mags_accessions.txt") as f:
    keep = set(line.strip() for line in f)
# prune
t.prune(keep, preserve_branch_length=True)

# write pruned tree
t.write(outfile="/home/roman.sanz/data/phylomodel/Clara_MAGs.tree",parser=0)
#format=0 writes internal nodel w/o node label, 
#prevents iToL from going crazy
