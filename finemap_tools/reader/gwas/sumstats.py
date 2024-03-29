from finemap_tools.reader import tabix_reader
import logging
import time
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
import pandas as pd
import numpy as np
from pathlib import Path
from finemap_tools.others.polyfun.polyfun_utils import set_snpid_index
from finemap_tools.utils import add_ID
from finemap_tools.snpfilter import filter_pipline
from scipy import stats


formatConvert2polyfun = {
    "pheweb": {"chrom": "CHR", "pos": "BP", "ref": "A2", "alt": "A1", "pval": "P"}
}


def formatConvert2polyfun_func(df, format):
    if format not in formatConvert2polyfun.keys():
        raise ValueError(
            f"format {format} is not supported, please check the formatConvert2polyfun"
        )

    df.rename(columns=formatConvert2polyfun[format], inplace=True)


def load_sumstats(
    sumstats_file,
    chr_num,
    start,
    end,
    allow_swapped_indel_alleles=False,
    sumstats_format=None,
):
    """
    sumstats_file should have these columns: SNP, CHR, BP, A1, A2, Z    (and optionally: P, SNPVAR)

    """
    # read sumstats and filter to target chromosome only
    logging.info("Loading sumstats file...")
    t0 = time.time()

    if sumstats_file.endswith(".gz") and Path(sumstats_file + ".tbi").exists():

        logging.info(
            f"loading with normal tabix reader to loading sumstats file with {sumstats_file}......"
        )
        df_sumstats = tabix_reader(sumstats_file, region=f"{chr_num}:{start}-{end}")

        logging.info(
            f"the header of loaded file is {df_sumstats.columns}, if it is not the correct header, please check the file or with correct format specified."
        )
        if all([isinstance(i, str) for i in df_sumstats.iloc[0].values]):
            logging.info(
                f"this tabix with header info in simply load and the first row is {df_sumstats.iloc[:1].values}"
            )
        if df_sumstats.shape[0] == 0:
            raise IOError(
                f"sumstats file does not include any SNPs in chromosome {chr_num} from {start} to {end}"
            )
        if sumstats_format is not None:
            formatConvert2polyfun_func(df_sumstats, sumstats_format)

    else:
        try:
            df_sumstats = pd.read_parquet(sumstats_file)
        except (ArrowIOError, ArrowInvalid):
            df_sumstats = pd.read_table(sumstats_file, sep="\s+")

        df_sumstats.rename(columns=formatConvert2polyfun[sumstats_format], inplace=True)
        if not np.any(df_sumstats["CHR"] == chr_num):
            raise IOError(
                "sumstats file does not include any SNPs in chromosome %s" % (chr_num)
            )
        formatConvert2polyfun_func(df_sumstats, sumstats_format)

        if np.any(df_sumstats["CHR"] != chr_num):
            df_sumstats = df_sumstats.query("CHR==%s" % (chr_num)).copy()

    if "SNP" not in df_sumstats.columns:
        df_sumstats["SNP"] = add_ID(df_sumstats, ["CHR", "BP", "A1", "A2"])

    df_sumstats = set_snpid_index(
        df_sumstats, allow_swapped_indel_alleles=allow_swapped_indel_alleles
    )

    if "P" not in df_sumstats.columns:
        df_sumstats["P"] = stats.chi2(1).sf(df_sumstats["Z"] ** 2)
    logging.info(
        "Loaded sumstats for %d SNPs in %0.2f seconds"
        % (df_sumstats.shape[0], time.time() - t0)
    )
    ## filter pipline provided by finemap_tools (if not installed will pass)

    logging.info(
        f"filtering SNP by finemap_tools with {df_sumstats.shape[0]} SNP at begining........"
    )
    df_sumstats = filter_pipline(sumstats=df_sumstats, id_col="SNP")

    logging.info(f"after filtering, left {df_sumstats.shape[0]} SNP")

    return df_sumstats
