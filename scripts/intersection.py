#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@Description:       :
@Date     :2024/03/05 15:19:30
@Author      :Tingfeng Xu
@version      :1.0
'''

"""
- eQTL:
SNP BETA SE GENE

- GWAS 
SNP BETA SE 

- pvar:
ID 

- gene expression martix 

1. check the used snps in eQTL and GWAS and pvar or your LD file

2. extract the used snps from all data files

"""




def parse_eQTL(file, ID):
    pass 

def parse_GWAS():
    pass 

def parse_pvar():
    pass 



import argparse
import sys
import warnings
import textwrap
from signal import SIG_DFL, SIGPIPE, signal


warnings.filterwarnings("ignore")
signal(
    SIGPIPE, SIG_DFL
)  # prevent IOError: [Errno 32] Broken pipe. If pipe closed by 'head'.


def getParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
        %prog Not work for big data! carefully when use large data 
        @Author: xutingfeng@big.ac.cn

        GWAS SSF: https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1
        GWAS Standard

        Version: 1.0
    
        """
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="input file",
        required=True, 
        default=[],
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="output file",
        required=True,
    )
    parser.add_argument(
        "-c", "--col", dest="col", help="column", required=True, default=[], nargs="+"
    )
    

    return parser
import pandas as pd 
if __name__ == "__main__":
    parser = getParser()
    args = parser.parse_args()
    
    inputFiles = args.input
    outputFile = args.output
    col_list = [eval(i) for i in args.col]

    set_list=[]
    for file, col in zip(inputFiles, col_list):
        # print(file, col)
        data = pd.read_csv(file, sep="\t", usecols=[col-1])

        data_set = set(data.iloc[:, 0].tolist())
        set_list.append(data_set)
        print(f"{file} used col is {data.columns[0]}: have {len(data_set)} rows")
    from functools import reduce
    res = reduce(lambda x, y: x.intersection(y), set_list)
    print(f"after merge, the total number of rows is {len(res)}")
    pd.DataFrame(res).to_csv(outputFile, index=False, header=False, sep="\t")


    #step1 find the used snps in eQTL and GWAS and pvar or your LD file
