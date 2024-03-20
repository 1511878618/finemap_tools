#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@Description:       :
@Date     :2024/03/19 16:42:43
@Author      :Tingfeng Xu
@version      :1.0
"""
"""
Source Code From: https://github.com/omerwe/polyfun
Modified By Tingfeng Xu
"""
from pathlib import Path
import argparse
import numpy as np
import textwrap
np.set_printoptions(precision=4, linewidth=200)
import pandas as pd

pd.set_option("display.width", 200)
import os
import logging

from finemap_tools.others.polyfun.polyfun_utils import (
    configure_logger,
    check_package_versions,
)
from finemap_tools.finemap import SUSIE_Wrapper, FINEMAP_Wrapper, coloc_Wrapper

# from polyfun import configure_logger, check_package_versions


def splash_screen():
    print("*********************************************************************")
    print("* Fine-mapping Wrapper")
    print("* Version 1.0.0")
    print("* (C) 2019-2022 Omer Weissbrod")
    print("*********************************************************************")
    print()


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

            
            """
        ),
    )

    # general parameters
    parser.add_argument(
        "--method",
        required=True,
        help="Fine-mapping method (currently susie and finemap are supported)",
        choices=["susie", "finemap", "GIFT", "coloc"],
    )
    parser.add_argument("--sumstats", required=True, help="Name of sumstats file")

    parser.add_argument("--chr", required=True, type=int, help="Target chromosome")
    parser.add_argument(
        "--start",
        required=True,
        type=int,
        help="First base-pair in the region to finemap",
    )
    parser.add_argument(
        "--end", required=True, type=int, help="Last base-pair in the region to finemap"
    )
    parser.add_argument("--n", required=True, type=int, help="Sample size")
    parser.add_argument(
        "--geno",
        default=None,
        help="Genotypes file (plink1 format, plink2 format(pgen) or bgen format)",
    )
    parser.add_argument(
        "--ld", default=None, help="prefix or fill name of an LD matrix file"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="name of the output file; while for coloc will work as output folder name",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="If specified, show verbose output",
    )
    parser.add_argument(
        "--debug-dir",
        default=None,
        help="If specified, this is a path of a directory that will include files for debugging problems",
    )
    parser.add_argument(
        "--sample-file",
        default=None,
        help="BGEN files must be used together with a sample file",
    )
    parser.add_argument(
        "--incl-samples",
        default=None,
        help="A single-column text file specifying the ids of individuals to include in fine-mapping",
    )

    # fine-mapping parameters
    parser.add_argument(
        "--max-num-causal", required=True, type=int, help="Number of causal SNPs"
    )
    parser.add_argument(
        "--non-funct",
        action="store_true",
        default=False,
        help="Perform non-functionally informed fine-mapping",
    )
    parser.add_argument(
        "--allow-missing",
        default=False,
        action="store_true",
        help="If specified, SNPs with sumstats that are not \
                            found in the LD panel will be omitted. This is not recommended, because the omitted SNPs may be causal,\
                            which could lead to false positive results",
    )
    parser.add_argument(
        "--allow-swapped-indel-alleles",
        default=False,
        action="store_true",
        help="If specified, indels whose alleles are swapped between the sumstats and LD matrix \
                            are kept for fine-mapping. The default behavior considers indels at the same position \
                            with swapped alleles to be different variants, and thus removes them. Use with caution. \
                            This is intended for use only when you are confident that the indels are identical, \
                            e.g. when using insample LD",
    )
    parser.add_argument(
        "--no-sort-pip",
        default=False,
        action="store_true",
        help="Do **not** sort results by PIP. Recommended for use with --susie-outfile",
    )

    # LDstore related parameters
    parser.add_argument(
        "--ldstore2", default=None, help="Path to an LDstore 2.0 executable file"
    )
    parser.add_argument(
        "--finemap-exe", default=None, help="Path to FINEMAP v1.4 executable file"
    )
    parser.add_argument(
        "--memory",
        type=int,
        default=1,
        help="Maximum amount of memory in GB to allocate to LDStore",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="The number of CPU cores LDstore will use (if not specified, LDstore will use the max number of CPU cores available",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="If specified, this is a path of a directory that will cache LD matrices that have already been computed",
    )
    parser.add_argument(
        "--cache-format",
        default=None,
        help="Format of the LDstore cache files (default: bcor); npz will with a metafile containing snpname named like .npz.gz should have columns: A1 (effect allele) A2 BP, SNP, chromosome",
        choices=["bcor", "npz", None],
    )
    # FINEMAP-specific parameters
    parser.add_argument(
        "--finemap-dir",
        default=None,
        help="If specified, the FINEMAP files will be saved to this directory",
    )

    # SuSiE-specific parameters
    parser.add_argument(
        "--susie-outfile",
        default=None,
        help="If specified, the SuSiE object will be saved to an output file (also see --no-sort-pip for help merging with main output file)",
    )
    parser.add_argument(
        "--susie-resvar",
        default=None,
        type=float,
        help="If specified, SuSiE will use this value of the residual variance",
    )
    parser.add_argument(
        "--susie-resvar-init",
        default=None,
        type=float,
        help="If specified, SuSiE will use this initial value of the residual variance",
    )
    parser.add_argument(
        "--susie-resvar-hess",
        default=False,
        action="store_true",
        help="If specified, SuSiE will specify the residual variance using the HESS estimate",
    )
    parser.add_argument(
        "--susie-max-iter",
        default=100,
        type=int,
        help="SuSiE argument max_iter which controls the max number of IBSS iterations to perform (default: 100)",
    )
    parser.add_argument(
        "--hess",
        action="store_true",
        default=False,
        help="If specified, estimate causal effect variance via HESS",
    )
    parser.add_argument(
        "--hess-iter",
        type=int,
        default=100,
        help="Average HESS over this number of iterations (default: 100)",
    )
    parser.add_argument(
        "--hess-min-h2",
        type=float,
        default=None,
        help="When estimating causal effect variance via HESS, exclude SNPs that tag less than this amount of heritability (default: None)",
    )
    # coloc
    ## fot qtl file
    parser.add_argument("--qtl", default=None, help="qtl file path, like GTExV8")
    parser.add_argument("--qtl-format", default="GTExV8", help="qtl format")
    parser.add_argument("--sumstats-format", default="format", help="sumstats format")
    parser.add_argument(
        "--gene_id", default="gene_id", help="gene_id col of qtl file to split"
    )
    # parser.add_argument("--n", default=None, help="number of samples in sumstats file")
    parser.add_argument(
        "--n2", default=None, type=int, help="number of samples in qtl file"
    )
    parser.add_argument(
        "--sdY1",
        default=None,
        type=int,
        help="sdY1, when --type1 is quant needed, while for cc, it is not needed",
    )
    # parser.add_argument("--sdY2", default=None, help="sdY2, when --type2 is quant needed, while for cc, it is not needed")
    parser.add_argument("--type1", default="quant", help="type1, quant or cc")
    # parser.add_argument("--type2", default="quant", help="type2, quant or cc")

    return parser.parse_args()


if __name__ == "__main__":

    args = args_parse()

    # check package versions
    check_package_versions()

    # show splash screen
    splash_screen()

    # check that the output directory exists
    if len(os.path.dirname(args.out)) > 0 and not os.path.exists(
        os.path.dirname(args.out)
    ):
        # raise ValueError(
        #     "output directory %s doesn't exist" % (os.path.dirname(args.out))
        # )

        Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    if (
        args.susie_outfile is not None
        and len(os.path.dirname(args.susie_outfile)) > 0
        and not os.path.exists(os.path.dirname(args.susie_outfile))
    ):
        # raise ValueError(
        #     "output directory %s doesn't exist" % (os.path.dirname(args.susie_outfile))
        # )
        Path(args.susie_outfile).parent.mkdir(parents=True, exist_ok=True)

    # configure logger
    configure_logger(args.out)

    # check params
    if args.max_num_causal == 1:
        if args.geno is not None or args.ld is not None:
            raise ValueError(
                "When max_num_causal=1 fine-mapping please omit the flags --geno and --ld (we cannot use LD information in that setting)"
            )
    else:
        if args.geno is None:
            if args.ld is None:
                raise ValueError("must specify either --geno or --ld")
            if args.ldstore2 is not None:
                raise ValueError("cannot specify both --ld and --ldstore2")
        if args.geno is not None:
            if args.ld is not None:
                raise ValueError("cannot specify both --geno and --ld")
            if args.geno.endswith(".bgen") and args.ldstore2 is None:
                raise ValueError(
                    "You must specify --ldstore2 when --geno that points to a bgen file"
                )

    if args.susie_outfile is not None and not args.no_sort_pip:
        logging.warning(
            "--susie-outfile was set but not --no-sort-pip. This will make it difficult to assign SNP names to the SuSiE R object"
        )

    # Create a fine-mapping class member
    if args.method == "susie":
        if args.finemap_dir is not None:
            raise ValueError("--finemap-dir cannot be specified with susie method")
        finemap_obj = SUSIE_Wrapper(
            genotypes_file=args.geno,
            sumstats_file=args.sumstats,
            n=args.n,
            chr_num=args.chr,
            start=args.start,
            end=args.end,
            sample_file=args.sample_file,
            incl_samples=args.incl_samples,
            ldstore_exe=args.ldstore2,
            n_threads=args.threads,
            cache_dir=args.cache_dir,
            cache_format=args.cache_format,
            memory=args.memory,
            allow_swapped_indel_alleles=args.allow_swapped_indel_alleles,
        )
    elif args.method == "finemap":
        if args.susie_outfile is not None:
            raise ValueError("--susie-outfile cannot be specified with finemap method")
        if args.finemap_exe is None:
            raise ValueError("need to specify --finemap-exe")
        if args.hess:
            raise ValueError("FINEMAP cannot be used with --hess")
        finemap_obj = FINEMAP_Wrapper(
            genotypes_file=args.geno,
            sumstats_file=args.sumstats,
            n=args.n,
            chr_num=args.chr,
            start=args.start,
            end=args.end,
            sample_file=args.sample_file,
            incl_samples=args.incl_samples,
            ldstore_exe=args.ldstore2,
            finemap_exe=args.finemap_exe,
            n_threads=args.threads,
            cache_dir=args.cache_dir,
            cache_format=args.cache_format,
            memory=args.memory,
            allow_swapped_indel_alleles=args.allow_swapped_indel_alleles,
        )
    elif args.method == "GIFT":
        # finemap_obj =

        raise NotImplementedError
    elif args.method == "coloc":

        finemap_obj = coloc_Wrapper(
            sumstats_file1=args.sumstats,
            sumstats_file1_format=args.sumstats_format,
            sumstats_file2=args.qtl,
            sumstats_file2_format=args.qtl_format,
            n1=args.n,
            n2=args.n2,
            type1=args.type1,
            sdY1=args.sdY1,
            chr_num=args.chr,
            start=args.start,
            end=args.end,
        )

    else:
        raise ValueError("unknown method specified in --method")

    # run fine-mapping
    if args.method == "coloc":
        finemap_obj.finemap(
            out_dir=args.out,
            gene_col=args.gene_id,
        )

    else:
        df_finemap = finemap_obj.finemap(
            locus_start=args.start,
            locus_end=args.end,
            num_causal_snps=args.max_num_causal,
            use_prior_causal_prob=not args.non_funct,
            prior_var=None,
            hess=args.hess,
            hess_iter=args.hess_iter,
            hess_min_h2=args.hess_min_h2,
            verbose=args.verbose,
            ld_file=args.ld,
            debug_dir=args.debug_dir,
            allow_missing=args.allow_missing,
            susie_outfile=args.susie_outfile,
            finemap_dir=args.finemap_dir,
            residual_var=args.susie_resvar,
            residual_var_init=args.susie_resvar_init,
            hess_resvar=args.susie_resvar_hess,
            susie_max_iter=args.susie_max_iter,
        )
        logging.info("Writing fine-mapping results to %s" % (args.out))
        if not args.no_sort_pip:
            df_finemap.sort_values("PIP", ascending=False, inplace=True)
        if args.out.endswith(".parquet"):
            df_finemap.to_parquet(args.out, index=False)
        else:
            df_finemap.to_csv(args.out, sep="\t", index=False, float_format="%0.5e")
