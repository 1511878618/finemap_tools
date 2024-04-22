#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@Description:       :
@Date     :2024/03/04 11:36:41
@Author      :Tingfeng Xu
@version      :1.0
"""


import argparse

import textwrap
import gwaslab as gl
import matplotlib.pyplot as plt
from finemap_tools.utils import add_ID
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s\n%(message)s",
    stream=sys.stdout,
)

from finemap_tools.utils import add_ID


def args_parse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
            %prog is ...
            @Author: xutingfeng@big.ac.cn
            Version: 1.0

            """
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="input gwas file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="output file suffix, output format is gwaslab",
    )
    parser.add_argument(
        "--format", default="regenie", required=True, help="input file format"
    )
    # parser.add_argument("--build", default="38", help="genome build version, only supported for 38")

    return parser.parse_args()


def main():
    args = args_parse()
    args = vars(args)

    sumstats = gl.Sumstats(args["input"], fmt=args["format"], build="38")
    if "SNPID" not in sumstats.data.columns:
        sumstats.data["SNPID"] = add_ID(
            sumstats.data, col_list=["CHR", "POS", "NEA", "EA"]
        )

        ## harmonize
        sumstats.harmonize(
            basic_check=True,
            n_cores=6,
            ref_seq=gl.get_path("ucsc_genome_hg38"),
            ref_infer=gl.get_path("1kg_eur_hg38"),
            ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg38_auto"),
            ref_alt_freq="AF",
            removedup_args=dict(remove=True, remove_dup=True),
        )
        from pathlib import Path

        Path(args["output"]).parent.mkdir(parents=True, exist_ok=True)
        sumstats.to_format(args["output"], bgzip=True, md5sum=True)


if __name__ == "__main__":
    main()
