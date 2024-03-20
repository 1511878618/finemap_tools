from .base import Fine_Mapping
from finemap_tools.utils import uri_validator, run_executable
import numpy as np


import pandas as pd

from finemap_tools.reader.ld import download_ld_file

import os
import time
import logging

import tempfile
from finemap_tools.others.polyfun.ldstore.bcor import bcor
from finemap_tools.reader.ld import download_ld_file, read_ld_from_file, get_bcor_meta


class FINEMAP_Wrapper(Fine_Mapping):

    def __init__(
        self,
        genotypes_file,
        sumstats_file,
        n,
        chr_num,
        start,
        end,
        finemap_exe,
        ldstore_exe,
        sample_file=None,
        incl_samples=None,
        cache_dir=None,
        cache_format=None,
        n_threads=None,
        memory=None,
        allow_swapped_indel_alleles=False,
    ):

        super(FINEMAP_Wrapper, self).__init__(
            genotypes_file,
            sumstats_file,
            n,
            chr_num=chr_num,
            start=start,
            end=end,
            ldstore_exe=ldstore_exe,
            sample_file=sample_file,
            incl_samples=incl_samples,
            cache_dir=cache_dir,
            cache_format=cache_format,
            n_threads=n_threads,
            memory=memory,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
        )
        self.finemap_exe = finemap_exe

    def finemap(
        self,
        locus_start,
        locus_end,
        num_causal_snps,
        use_prior_causal_prob=True,
        prior_var=None,
        residual_var=None,
        hess=False,
        hess_iter=100,
        hess_min_h2=None,
        susie_max_iter=100,
        verbose=False,
        ld_file=None,
        debug_dir=None,
        allow_missing=False,
        susie_outfile=None,
        residual_var_init=None,
        hess_resvar=False,
        finemap_dir=None,
    ):

        # check params
        if use_prior_causal_prob and "SNPVAR" not in self.df_sumstats.columns:
            raise ValueError("SNPVAR column not found in sumstats file")
        if hess:
            raise ValueError(
                "FINEMAP cannot be used with a HESS-based variance estimator"
            )
        if residual_var is not None:
            raise ValueError("cannot specify residual_var for FINEMAP")
        if debug_dir is not None:
            raise NotImplementedError("FINEMAP object does not support --debug-dir")
        if hess_resvar:
            raise NotImplementedError(
                "FINEMAP object does not support --susie-resvar-hess"
            )
        if residual_var_init is not None:
            raise NotImplementedError(
                "FINEMAP object does not support --susie-resvar-init"
            )
        # if allow_missing:
        # raise ValueError('FINEMAP object does not support --allow-missing')

        # download LD file if it's a url
        if uri_validator(ld_file):
            ld_file = download_ld_file(ld_file)

        # create prefix of output files
        if finemap_dir is None:
            finemap_dir = tempfile.mkdtemp()
        else:
            os.makedirs(finemap_dir, exist_ok=True)
            logging.info("Saving FINEMAP files to directory: %s" % (finemap_dir))
        assert os.path.isdir(finemap_dir)
        finemap_output_prefix = os.path.join(finemap_dir, "finemap")

        # set locus
        self.set_locus(locus_start, locus_end)

        # find or create a suitable ld_file
        if num_causal_snps == 1:
            if ld_file is not None:
                raise ValueError(
                    "cannot specify an ld file when assuming a single causal SNP per locus"
                )
            ld_file = finemap_output_prefix + ".ld"
            np.savetxt(
                ld_file,
                np.eye(self.df_sumstats_locus.shape[0], dtype=np.int64),
                fmt="%s",
            )
        else:
            if ld_file is None:
                ld_data = self.get_ld_data(
                    locus_start, locus_end, need_bcor=True, verbose=verbose
                )
                if isinstance(ld_data, str):
                    ld_file = ld_data
                    assert ld_file.endswith(".bcor")
                    assert os.path.exists(ld_file)
                elif isinstance(ld_data, tuple):
                    assert len(ld_data) == 2
                    ld_arr, df_ld_snps = ld_data[0], ld_data[1]
                    self.sync_ld_sumstats(
                        ld_arr, df_ld_snps, allow_missing=allow_missing
                    )
                    del ld_arr, df_ld_snps
                    ld_file = finemap_output_prefix + ".ld"
                    np.savetxt(ld_file, self.df_ld.values, fmt="%0.5f")
            elif ld_file.endswith(".bcor"):
                pass
            elif ld_file.endswith(".npz") or os.path.exists(ld_file + ".npz"):
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
                del ld_arr, df_ld_snps
                ld_file = finemap_output_prefix + ".ld"
                np.savetxt(ld_file, self.df_ld.values, fmt="%0.5f")
            else:
                raise ValueError("unknown LD file format for file: %s" % (ld_file))

        # define file names
        master_file = finemap_output_prefix + ".master"
        snp_filename = finemap_output_prefix + ".snp"
        config_filename = finemap_output_prefix + ".config"
        cred_filename = finemap_output_prefix + ".cred"
        log_filename = finemap_output_prefix + ".log"
        z_filename = finemap_output_prefix + ".z"

        # flip some of the alleles
        if num_causal_snps == 1:
            is_flipped = np.zeros(self.df_sumstats_locus.shape[0], dtype=bool)
        else:
            if ld_file.endswith(".bcor"):
                bcor_obj = bcor(ld_file)
                df_ld_snps = get_bcor_meta(bcor_obj)
                self.sync_ld_sumstats(None, df_ld_snps, allow_missing=allow_missing)
                del df_ld_snps
            assert np.all(self.df_ld_snps["BP"] == self.df_sumstats_locus["BP"])
            is_flipped = self.df_ld_snps["A1"] == self.df_sumstats_locus["A2"]
            is_not_flipped = self.df_ld_snps["A1"] == self.df_sumstats_locus["A1"]
            assert np.all(is_flipped | is_not_flipped)
            if np.any(is_flipped):
                logging.info(
                    "Flipping the effect-sign of %d SNPs that are flipped compared to the LD panel"
                    % (is_flipped.sum())
                )

        # create df_z and save it to disk
        df_z = self.df_sumstats_locus[["SNP", "CHR", "BP", "A1", "A2", "Z"]].copy()
        df_z.loc[is_flipped, "A1"] = self.df_sumstats_locus.loc[is_flipped, "A2"]
        df_z.loc[is_flipped, "A2"] = self.df_sumstats_locus.loc[is_flipped, "A1"]
        df_z.loc[is_flipped, "Z"] *= -1

        df_z.rename(
            columns={
                "SNP": "rsid",
                "CHR": "chromosome",
                "BP": "position",
                "A1": "allele1",
                "A2": "allele2",
                "Z": "beta",
            },
            inplace=True,
            errors="ignore",
        )
        df_z["se"] = 1
        df_z["maf"] = 0.05
        df_z = df_z[
            [
                "rsid",
                "chromosome",
                "position",
                "allele1",
                "allele2",
                "maf",
                "beta",
                "se",
            ]
        ]
        if use_prior_causal_prob:
            df_z["prob"] = (
                self.df_sumstats_locus["SNPVAR"]
                / self.df_sumstats_locus["SNPVAR"].sum()
            )
        df_z.to_csv(z_filename, header=True, index=False, sep=" ", float_format="%0.5f")

        # create the master file
        df_master = pd.DataFrame()
        df_master["z"] = [z_filename]
        df_master["snp"] = [snp_filename]
        df_master["config"] = [config_filename]
        df_master["cred"] = [cred_filename]
        df_master["log"] = [log_filename]
        df_master["n_samples"] = [self.n]
        if ld_file.endswith(".bcor"):
            df_master["bcor"] = [ld_file]
        elif ld_file.endswith(".ld"):
            df_master["ld"] = [ld_file]
        else:
            raise ValueError("Illegal LD file format")
        df_master.to_csv(master_file, sep=";", header=True, index=False)

        # prepare the FINEMAP command
        finemap_cmd = [self.finemap_exe]
        finemap_cmd += ["--in-files", master_file, "--sss"]
        finemap_cmd += ["--force-n-samples"]
        finemap_cmd += ["--log", log_filename]
        finemap_cmd += ["--n-causal-snps", str(num_causal_snps)]
        finemap_cmd += ["--std-effects"]
        ###finemap_cmd += ['--flip-beta']
        if self.n_threads is not None:
            finemap_cmd += ["--n-threads", str(self.n_threads)]
        if prior_var is not None:
            finemap_cmd += ["--prior-std", str(np.sqrt(prior_var))]
        if use_prior_causal_prob:
            finemap_cmd += ["--prior-snps"]

        # run FINEMAP
        t0 = time.time()
        m = self.df_sumstats_locus.shape[0]
        logging.info(
            "Starting %s FINEMAP fine-mapping for chromosome %d BP %d-%d (%d SNPs)"
            % (
                (
                    "functionally-informed"
                    if use_prior_causal_prob
                    else "non-functionally informed"
                ),
                self.chr,
                locus_start,
                locus_end,
                self.df_sumstats_locus.shape[0],
            )
        )
        run_executable(
            finemap_cmd,
            "FINEMAP",
            measure_time=True,
            show_output=verbose,
            show_command=verbose,
        )
        if not os.path.exists(log_filename + "_sss"):
            raise IOError("FINEMAP output files not found")

        # load log file
        found_post_csnps = False
        with open(log_filename + "_sss") as f:
            for line in f:
                if line.startswith("- Post-expected # of causal SNPs"):
                    post_mean_num_csnps = float(line[line.index(": ") + 2 : -1])
                    found_post_csnps = True
                    break
        if not found_post_csnps:
            raise IOError("corrupt log file found: %s" % (log_filename + "_sss"))

        # load results
        df_finemap = pd.read_table(
            snp_filename,
            sep=" ",
            usecols=[
                "rsid",
                "chromosome",
                "position",
                "allele1",
                "allele2",
                "prob",
                "mean",
                "sd",
            ],
        )
        df_finemap.rename(
            columns={
                "rsid": "SNP",
                "position": "BP",
                "chromosome": "CHR",
                "prob": "PIP",
                "mean": "BETA_MEAN",
                "sd": "BETA_SD",
                "allele1": "A1",
                "allele2": "A2",
            },
            inplace=True,
            errors="raise",
        )

        # read log10bf
        log10bf = None
        with open(log_filename + "_sss", "r") as f:
            for line in f:
                if line.startswith("- Log10-BF"):
                    log10bf = float(line.split()[-1])
                    break
        if log10bf is None:
            raise ValueError("FINEMP did not report Log10-BF")

        # add distance from center
        start = df_finemap["BP"].min()
        end = df_finemap["BP"].max()
        middle = (start + end) // 2
        df_finemap["DISTANCE_FROM_CENTER"] = np.abs(df_finemap["BP"] - middle)

        # add causal set info
        df_finemap["CREDIBLE_SET"] = 0
        cred_file = None
        for m in range(num_causal_snps, 0, -1):
            if os.path.exists(cred_filename + str(m)):
                cred_file = cred_filename + str(m)
                break
        if cred_file is None:
            raise IOError("cred file not found")
        df_cred = pd.read_table(
            cred_file, sep=" ", usecols=(lambda c: c.startswith("cred")), comment="#"
        )
        df_finemap.set_index("SNP", inplace=True, drop=False)
        for c_i, c in enumerate(df_cred.columns):
            df_finemap.loc[df_cred[c].dropna(), "CREDIBLE_SET"] = c_i + 1
        df_finemap.reset_index(inplace=True, drop=True)

        finemap_time = time.time() - t0
        logging.info("Done in %0.2f seconds" % (finemap_time))

        return df_finemap
