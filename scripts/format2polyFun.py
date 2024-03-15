#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@Description:       :
@Date     :2024/03/04 11:36:41
@Author      :Tingfeng Xu
@version      :1.0
"""

try:
    import pandas as pd
except:
    import subprocess
    import sys

    print("installing pandas")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "pandas"], check=True)
        print("pandas模块安装完成。")
        import pandas as pd
    except subprocess.CalledProcessError:
        print("安装pandas模块时出错。请手动安装pandas。")
        sys.exit(1)

import argparse
import os
import textwrap
import warnings
import logging
import sys
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s\n%(message)s', stream=sys.stdout)

from finemap_tools.utils import add_ID

def args_parse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
            %prog is ...
            @Author: xutingfeng@big.ac.cn
            Version: 1.0
            Turn the sumstats to polyFun format
            Only support for all string cols name or all int cols index as input 
            Example code:
                format2polyFun.py -i DNAJC16.tsv -o test.tsv.gz --chr chrom --pos pos --A1 alt --A2 ref --beta BETA --se SE  
            """
        ),
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="input sumstats, will use pandas to read the file")
    parser.add_argument("-s", "--sep", type=str, required=False, default="\t", help="sep of the input file, default is '\\t' ")
    parser.add_argument("-o", "--output", type=str, required=True, help="output folder, if with .gz will be compressed")
    parser.add_argument("--snp", type=str, required=False, help="snp col, if not given, will use the chrom:pos:other_allele:effect_allele as snp col")
    parser.add_argument("--chr", type=str, required=True, help="chr col")
    parser.add_argument("--pos", type=str, required=True, help="pos col")
    parser.add_argument("--A1", type=str, required=True, help="A1 is the effect allele (i.e. the sign of Z should match the A1 allele)")
    parser.add_argument("--A2", type=str, required=True, help="A2 is the other allele (i.e. the reference allele)")
    parser.add_argument("--Z", type=str, required=False, help="Z is the signed Z-score for the effect allele, if not passed will use the effect size and standard error to calculate Z, so make sure passed with --beta and --se and the effect size is the effect size of A1")
    parser.add_argument("--beta", type=str, required=False, help="beta is the effect size of the effect allele")
    parser.add_argument("--se", type=str, required=False, help="se is the standard error of the effect size")
    parser.add_argument("--no-header", action="store_true", help="if passed, will not use the first row as header")
    # parser.add_argument("--parallel", action="store_true", help="if passed, will use parallel to split the file")
    return parser.parse_args()


def format2polyFun(data, args ):

    if not args['snp']:
        data['SNP'] = add_ID(data, [args['chr'], args['pos'], args['A2'], args['A1']])
    else:
        data.rename(columns={args['snp']: "SNP"}, inplace=True)
    
    if not args['Z']:
        if not args['beta'] and not args['se']:
            raise ValueError("if --Z not passed, --beta and --se should be passed")
        data['Z'] = data[args['beta']] / data[args['se']]
    else:
        data.rename(columns={args['Z']: "Z"}, inplace=True)

    data = data.rename(columns={args['chr']: "CHR", args['pos']: "BP", args['A1']: "A1", args['A2']: "A2"})    
    return data[['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z']]

def map_header(header, col_idx):
    return header[col_idx-1]

def is_int(s):
    try:
        int(s)
        return True
    except:
        return False

def main():
    args = args_parse()
    args = vars(args)

    if not args['snp']:
        logging.info("snp col not passed, will use the chrom:pos:other_allele:effect_allele as snp col, make sure the rsid of bgen file is consistent with the chrom:pos:other_allele:effect_allele")

    if not args['Z']:
        if not args['beta'] and not args['se']:
            raise ValueError("if --Z not passed, --beta and --se should be passed")
        logging.info(
            "calculating Z from beta and se, make sure the beta is consistent with the effect allele!"
        )
    default_polyFun_cols = ['snp', 'chr', 'pos', 'A1', 'A2', 'Z', "beta", "se"]
    use_cols = [
        int(args[k]) if is_int(args[k]) else args[k]
        for k in default_polyFun_cols
        if args[k]
    ]
    data = pd.read_csv(
        args["input"],
        sep=args["sep"],
        header=None if args["no_header"] else 0,
        usecols=use_cols,
    )
    for key in default_polyFun_cols:
        if is_int(args[key]): # 如果是数字，就是列索引,start from 1
            args[key] = map_header(data.columns, int(args[key]))
    logging.info(f"args:\n {args}")

    data = format2polyFun(data, args)

    data.to_csv(args['output'], sep="\t", index=False, compression="gzip" if args['output'].endswith(".gz") else None)
    logging.info(f"output to {args['output']} and done!")


if __name__ == "__main__":
    main()
