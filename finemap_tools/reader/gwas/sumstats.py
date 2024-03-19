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


def load_sumstats(sumstats_file, chr_num, allow_swapped_indel_alleles=False):
    """
    sumstats_file should have these columns: SNP, CHR, BP, A1, A2, Z    (and optionally: P, SNPVAR)

    """
    # read sumstats and filter to target chromosome only
    logging.info("Loading sumstats file...")
    t0 = time.time()

    if sumstats_file.endswith(".gz") and Path(sumstats_file + ".tbi").exists():

        df_sumstats = tabix_reader(sumstats_file, region=f"{chr_num}")
        if all([isinstance(i, str) for i in df_sumstats.iloc[0].values]):
            logging.info(
                f"this tabix with header info in simply load and the first row is {df_sumstats.iloc[:1].values}"
            )

    else:
        try:
            df_sumstats = pd.read_parquet(sumstats_file)
        except (ArrowIOError, ArrowInvalid):
            df_sumstats = pd.read_table(sumstats_file, sep="\s+")
        if not np.any(df_sumstats["CHR"] == chr_num):
            raise IOError(
                "sumstats file does not include any SNPs in chromosome %s" % (chr_num)
            )
        if np.any(df_sumstats["CHR"] != chr_num):
            df_sumstats = df_sumstats.query("CHR==%s" % (chr_num)).copy()

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

    df_sumstats["added_id"] = add_ID(df_sumstats, ["CHR", "BP", "A1", "A2"])

    logging.info(
        f"filtering SNP by finemap_tools with {df_sumstats.shape[0]} SNP at begining........"
    )
    df_sumstats = filter_pipline(sumstats=df_sumstats, id_col="added_id")

    logging.info(f"after filtering, left {df_sumstats.shape[0]} SNP")
    df_sumstats = df_sumstats.drop(columns=["added_id"])
    return df_sumstats
