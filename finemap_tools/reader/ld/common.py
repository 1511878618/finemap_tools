import logging
import time
import os
import pandas as pd
from scipy import sparse
import numpy as np
from finemap_tools.others.polyfun.ldstore.bcor import bcor
import numpy as np


import pandas as pd


import os
import time
import logging

import tempfile
from finemap_tools.others.polyfun.polyfun_utils import (
    TqdmUpTo,
)
from finemap_tools.others.polyfun.ldstore.bcor import bcor
import scipy.sparse as sparse

# from polyfun import configure_logger, check_package_versions
import urllib.request


def load_ld_npz(ld_prefix):

    logging.info("Loading LD file %s" % (ld_prefix))
    t0 = time.time()

    # load SNPs info
    snps_filename_parquet = ld_prefix + ".parquet"
    snps_filename_gz = ld_prefix + ".gz"
    if os.path.exists(snps_filename_parquet):
        df_ld_snps = pd.read_parquet(snps_filename_parquet)
    elif os.path.exists(snps_filename_gz):
        df_ld_snps = pd.read_table(snps_filename_gz, sep="\s+")

        # df_ld_snps.rename(columns={'allele1':'A1', 'allele2':'A2', 'position':'BP', 'chromosome':'CHR', 'rsid':'SNP'}, inplace=True, errors='ignore')
    else:
        raise ValueError(
            "couldn't find SNPs file %s or %s"
            % (snps_filename_parquet, snps_filename_gz)
        )

    # load LD matrix
    R_filename = ld_prefix + ".npz"
    if not os.path.exists(R_filename):
        raise IOError("%s not found" % (R_filename))
    ld_arr = sparse.load_npz(R_filename).toarray()
    ld_arr = ld_arr + ld_arr.T
    assert np.allclose(np.diag(ld_arr), 1.0)
    # assert np.all(~np.isnan(ld_arr))

    # sanity checks
    assert ld_arr.shape[0] == ld_arr.shape[1]
    if ld_arr.shape[0] != df_ld_snps.shape[0]:
        raise ValueError("LD matrix has a different number of SNPs than the SNPs file")

    logging.info("Done in %0.2f seconds" % (time.time() - t0))
    return ld_arr, df_ld_snps


def get_bcor_meta(bcor_obj):
    df_ld_snps = bcor_obj.getMeta()
    df_ld_snps.rename(
        columns={
            "rsid": "SNP",
            "position": "BP",
            "chromosome": "CHR",
            "allele1": "A1",
            "allele2": "A2",
        },
        inplace=True,
        errors="raise",
    )
    ###df_ld_snps['CHR'] = df_ld_snps['CHR'].astype(np.int64)
    df_ld_snps["BP"] = df_ld_snps["BP"].astype(np.int64)
    return df_ld_snps


def load_ld_bcor(ld_prefix):
    bcor_file = ld_prefix + ".bcor"
    if not os.path.exists(bcor_file):
        raise IOError("%s not found" % (bcor_file))
    logging.info("Loading LD file %s" % (bcor_file))
    t0 = time.time()
    bcor_obj = bcor(bcor_file)
    df_ld_snps = get_bcor_meta(bcor_obj)
    ld_arr = bcor_obj.readCorr([])
    # assert np.all(~np.isnan(ld_arr))
    logging.info("Done in %0.2f seconds" % (time.time() - t0))
    return ld_arr, df_ld_snps


def read_ld_from_file(ld_file):
    # if ld_file is a prefix, make it into a full file name
    if not ld_file.endswith(".bcor") and not ld_file.endswith(".npz"):
        if os.path.exists(ld_file + ".npz"):
            ld_file = ld_file + ".npz"
        elif os.path.exists(ld_file + ".bcor"):
            ld_file = ld_file + ".bcor"
        else:
            raise IOError("No suitable LD file found")

    # read the LD file
    if ld_file.endswith(".bcor"):
        ld_arr, df_ld_snps = load_ld_bcor(ld_file[:-5])  # TODO: modify
    elif ld_file.endswith(".npz"):
        ld_arr, df_ld_snps = load_ld_npz(ld_file[:-4])  # TODO:modify
    else:
        raise ValueError("unknown LD format")
    # is_na_ld = np.all(
    #     ~np.isnan(ld_arr)
    # )  # only keep this and I suppose this check is no need, only thing could do is to avoid the nan in the ld_arr, and drop them all.
    logging.warning(
        f"there are {np.isnan(ld_arr).sum(axis=1).shape[0]} nan in R matrix in the ld_arr"
    )

    return ld_arr, df_ld_snps


def download_ld_file(url_prefix):
    temp_dir = tempfile.mkdtemp()
    filename_prefix = os.path.join(temp_dir, "ld")
    for suffix in ["npz", "gz"]:
        url = url_prefix + "." + suffix
        suffix_file = filename_prefix + "." + suffix
        with TqdmUpTo(
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            miniters=1,
            desc="downloading %s" % (url),
        ) as t:
            try:
                urllib.request.urlretrieve(
                    url, filename=suffix_file, reporthook=t.update_to
                )
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    raise ValueError("URL %s wasn't found" % (url))
                else:
                    raise

    return filename_prefix
