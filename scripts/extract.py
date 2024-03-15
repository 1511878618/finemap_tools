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
        %prog 
        @Author: xutingfeng@big.ac.cn

        GWAS SSF: https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1
        GWAS Standard

        Version: 1.0
    
        """
        ),
    )

    parser.add_argument('--extract', dest='extract', help='extract by given files, this file should only have 1 columns', default=None, required=True)
    parser.add_argument('-c', '--col', dest='col', help='column number, start from 1 and default is 1', default=1, required=False, type=int)
    parser.add_argument('-d', '--delimiter', dest='delimiter', help='delimiter, default is tab', default=None, required=False)
    parser.add_argument('--no-header', dest='no_header', help='no header', action='store_true', default=False, required=False)


    return parser
import pandas as pd 
if __name__ == "__main__":
    parser = getParser()
    args = parser.parse_args()

    delimiter = args.delimiter
    extract = args.extract
    col = args.col
    extract_list = pd.read_csv(extract, header=None, sep='\t', usecols=[0]).iloc[:, 0].tolist()

    out_delimiter = '\t' if delimiter is None else delimiter
    line_idx = 1 
    for line in sys.stdin:
        if line_idx == 1 and not args.no_header:
            sys.stdout.write(line)
            line_idx += 1
            continue

        line = line.strip()
        line_split = line.split(delimiter)        
        if line_split[col-1] in extract_list:
            sys.stdout.write(out_delimiter.join(line_split) + '\n')
        line_idx += 1

    sys.stdout.close()
    sys.stderr.flush()
    sys.stderr.close()
       



