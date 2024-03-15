gwasPath = "/home/xutingfeng/GIFT/data/GWAS/DNAJC16_gift.tsv" 
eQTLsLocation = "/home/xutingfeng/GIFT/data/eQTL/Kidney_Cortex/genes"
LDPath = "/home/xutingfeng/GIFT/data/pgen/DNAJC16_LD.unphased.vcor2"
snplistPath = "/home/xutingfeng/GIFT/data/pgen/DNAJC16_LD.unphased.vcor2.vars"
geneExpression = "/home/xutingfeng/GIFT/data/expression/gene_tpm_2017-06-05_v8_kidney_cortex.gct.gz"



from pathlib import Path
import numpy as np 
from functools import reduce

def load_eQTLs(eQTLsLocation, snplist, gene_list=None, SNP="ID", BETA="BETA", SE="SE" ):
    gene_path_dict = {i.stem: str(i) for i in Path(eQTLsLocation).glob("*") }
    print(f"found gene files: {len(gene_path_dict)}")

    if gene_list is not None:
        tmp_dict = {}
        for gene in gene_path_dict:
            if gene in gene_list:
                tmp_dict[gene] = gene_path_dict[gene]
        print(f"only {len(tmp_dict)} in passed gene list ")
        gene_path_dict = tmp_dict       

    gene_list = []
    pindex = []
    final = [] 
    snplist_df = pd.DataFrame({SNP: snplist})
    for gene, path in gene_path_dict.items():
        eQTLs = pd.read_csv(path, sep='\t').drop_duplicates(subset=[SNP])
        eQTLs[f'{gene}'] = eQTLs[BETA] / eQTLs[SE]
        
        # incase for one eQTLs may have duplicate SNPs as in different ProbeID
        eQTLs = snplist_df.merge(eQTLs, on=SNP, how='left')[[SNP, f"{gene}"]]
        final.append(eQTLs)

        gene_list.append(gene)
        pindex.append( eQTLs[f"{gene}"].notnull().sum())
    return reduce(lambda x, y: pd.merge(x, y, on=SNP, how='inner'), final).set_index(SNP), gene_list, pindex


import pandas as pd 
geneExpression_df = pd.read_csv(geneExpression, sep = '\t', skiprows=2 ).iloc[:, 2:].set_index('Description').T
# 去除重复的基因
geneExpression_df = geneExpression_df.loc[:, ~geneExpression_df.columns.duplicated()]

all_gene_list = geneExpression_df.columns.tolist()
all_gene_list

snplist = pd.read_csv(snplistPath, sep='\t', header=None).iloc[:, 0].tolist()
eQTLs, gene_list, pindex= load_eQTLs(eQTLsLocation,snplist, all_gene_list)

eQTLs

gwas = pd.read_csv(gwasPath, sep='\t')
gwas['TSTAT'] = gwas['BETA'] / gwas['SE']
snplist_df = pd.DataFrame({"ID":snplist})
gwas = snplist_df.merge(gwas, on="ID", how='left').set_index("ID")[['TSTAT']]
gwas

LD = pd.read_csv(LDPath, sep='\t', header=None)
LD 

R = geneExpression_df[gene_list].corr() 
R 

assert eQTLs.columns.tolist() == gene_list
assert eQTLs.index.tolist() == snplist
assert gwas.index.tolist() == snplist
assert R.index.tolist() == gene_list
assert R.columns.tolist() == gene_list

import h5py
with h5py.File('test.h5', 'w') as f:
    f.create_group("data")
    f['data'].create_dataset('eQTLs', data=eQTLs.values)
    f['data'].create_dataset('gwas', data=gwas.values)
    f['data'].create_dataset('LD', data=LD.values)
    f['data'].create_dataset('R', data=R.values)
    f['data'].create_dataset('pindex', data=pindex)
    f['data'].create_dataset('gene', data=[i.encode() for i in gene_list])

    # f['eQTLs'] = eQTLs.values # row is SNP, column is gene
    # f['eQTLsLD'] = LD.values # row is SNP, column is SNP
    # f['gawsLD'] = LD.values #  row is SNP, column is SNP
    # f['gwas'] = gwas.values # row is SNP, column is TSTAT
    # f['pindex'] = pindex # row is gene
