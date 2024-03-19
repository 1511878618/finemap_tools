from .base import Fine_Mapping
import numpy as np
import pandas as pd

import os
import time
import logging

import shutil
import glob
from importlib import reload
from packaging.version import Version
from finemap_tools.utils import uri_validator
from finemap_tools.reader.ld import download_ld_file, read_ld_from_file


class SUSIE_Wrapper(Fine_Mapping):

    def __init__(
        self,
        genotypes_file,
        sumstats_file,
        n,
        chr_num,
        ldstore_exe,
        sample_file=None,
        incl_samples=None,
        cache_dir=None,
        cache_format=None,
        n_threads=None,
        memory=None,
        allow_swapped_indel_alleles=False,
    ):

        super(SUSIE_Wrapper, self).__init__(
            genotypes_file,
            sumstats_file,
            n,
            chr_num,
            ldstore_exe=ldstore_exe,
            sample_file=sample_file,
            incl_samples=incl_samples,
            cache_dir=cache_dir,
            cache_format=cache_format,
            n_threads=n_threads,
            memory=memory,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
        )

        # load SuSiE R package
        import rpy2
        import rpy2.robjects.numpy2ri as numpy2ri
        import rpy2.robjects as ro

        ro.conversion.py2ri = numpy2ri
        numpy2ri.activate()
        from rpy2.robjects.packages import importr

        self.susieR = importr("susieR")
        self.R_null = ro.rinterface.NULL
        # self.RNULLType = rpy2.rinterface.RNULLType

    def finemap(
        self,
        locus_start,
        locus_end,
        num_causal_snps,
        use_prior_causal_prob=True,
        prior_var=None,
        residual_var=None,
        residual_var_init=None,
        hess_resvar=False,
        hess=False,
        hess_iter=100,
        hess_min_h2=None,
        susie_max_iter=100,
        verbose=False,
        ld_file=None,
        debug_dir=None,
        allow_missing=False,
        susie_outfile=None,
        finemap_dir=None,
    ):

        # check params
        if use_prior_causal_prob and "SNPVAR" not in self.df_sumstats.columns:
            raise ValueError("SNPVAR column not found in sumstats file")
        if hess_resvar:
            assert hess, "hess_resvar cannot be specified if hess is FALSE"

        # set locus
        self.set_locus(locus_start, locus_end)

        # download LD file if it's a url
        if uri_validator(ld_file):
            ld_file = download_ld_file(ld_file)
            delete_ld_files_on_exit = True
        else:
            delete_ld_files_on_exit = False

        # Load LD data into memory if num_causal_snps>1
        if num_causal_snps == 1:
            if hess:
                raise ValueError(
                    "Cannot use HESS-based variance estimator when assuming a single causal SNP per locus"
                )
            self.df_ld = pd.DataFrame(
                np.eye(self.df_sumstats_locus.shape[0]),
                index=self.df_sumstats_locus.index,
                columns=self.df_sumstats_locus,
            )
            self.df_ld_snps = self.df_sumstats_locus
        else:
            if ld_file is None:
                ld_arr, df_ld_snps = self.get_ld_data(
                    locus_start, locus_end, need_bcor=False, verbose=verbose
                )
            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)

            # assert np.all(~np.isnan(ld_arr)) # NOTE: this may often happned with some ld is NaN, I supposed to rm them instead of raise errors; And code in self.sync_ld_sumstats will do this, so i comments this code here

            self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
            del ld_arr
            del df_ld_snps

        # define prior causal probabilities
        if use_prior_causal_prob:
            prior_weights = self.df_sumstats_locus["SNPVAR"].copy().values
            prior_weights /= prior_weights.sum()
            assert np.isclose(prior_weights.sum(), 1)

        # flip effect sizes if needed
        assert np.all(self.df_ld_snps["BP"] == self.df_sumstats_locus["BP"])
        is_flipped = self.df_ld_snps["A1"] == self.df_sumstats_locus["A2"]
        is_not_flipped = self.df_ld_snps["A1"] == self.df_sumstats_locus["A1"]
        assert np.all(is_flipped | is_not_flipped)
        bhat = self.df_sumstats_locus["Z"].values.copy()
        if np.any(is_flipped):
            bhat[is_flipped.values] *= -1
            logging.info(
                "Flipping the effect-sign of %d SNPs that are flipped compared to the LD panel"
                % (is_flipped.sum())
            )

        # Use HESS to estimate causal effect sizes
        if hess:
            if prior_var is not None:
                raise ValueError("cannot specify both hess and a custom prior_var")
            if self.n < 20000:
                logging.warning(
                    "HESS method is intended for studies with large sample sizes (i.e. >20K)"
                )
            if hess_min_h2 is None:
                logging.warning(
                    "For best results, you should consider setting --hess-min-h2 to exclude SNPs with low heritability from the HESS estimation. You will need to experiment with your data to find a suitable heritability threshold. To start, try --hess-min-h2 1e-4"
                )
            else:
                logging.info(
                    "Excluding SNPs with heritability less than %0.4e from the HESS estimation"
                    % (hess_min_h2)
                )
            h2_hess = self.estimate_h2_hess_wrapper(
                min_h2=hess_min_h2, num_samples=hess_iter
            )
            logging.info(
                "Average local SNP heritability estimated by modified HESS over %d iterations: %0.4e"
                % (hess_iter, h2_hess)
            )
            if h2_hess > 10:
                logging.warning(
                    "The HESS estimator is unconstrained, and the estimate is an order of magnitude greater than the expected max of 1. Use with caution"
                )
            prior_var = h2_hess / num_causal_snps
            if prior_var <= 0:
                raise ValueError(
                    "HESS estimates that the locus causally explains zero heritability"
                )
            if prior_var >= 1:
                raise ValueError(
                    "HESS-estimated prior-var >1. The HESS estimator cannot be used in this locus."
                )
            logging.info(
                "HESS estimated causal effect size variance: %0.4e" % (prior_var)
            )

            if hess_resvar:
                residual_var = 1 - h2_hess
                logging.info(
                    "Residual variance using the HESS estimate: %0.4e" % (residual_var)
                )
                assert residual_var >= 0

        # rpy2 bug fix
        import rpy2.robjects.numpy2ri as numpy2ri

        reload(numpy2ri)
        numpy2ri.activate()

        # run SuSiE
        t0 = time.time()
        m = self.df_sumstats_locus.shape[0]
        logging.info(
            "Starting %s SuSiE fine-mapping for chromosome %d BP %d-%d (%d SNPs)"
            % (
                (
                    "functionally-informed"
                    if use_prior_causal_prob
                    else "non-functionally informed"
                ),
                self.chr,
                locus_start,
                locus_end,
                self.df_ld.shape[0],
            )
        )

        # save variables to debug dir if needed
        if debug_dir is not None:
            os.makedirs(debug_dir, exist_ok=True)
            logging.info("Saving debug info to: %s" % (debug_dir))
            self.df_sumstats_locus.to_csv(
                os.path.join(debug_dir, "df_sumstats_locus.txt"), index=False, sep="\t"
            )
            np.savetxt(os.path.join(debug_dir, "bhat.txt"), bhat)
            # np.savez_compressed(os.path.join(debug_dir, 'R.npz'), R=self.df_ld.values)
            np.savetxt(os.path.join(debug_dir, "n.txt"), [self.n])
            np.savetxt(os.path.join(debug_dir, "L.txt"), [num_causal_snps])
            np.savetxt(
                os.path.join(debug_dir, "residual_var.txt"),
                [np.nan] if (residual_var is None) else [residual_var],
            )
            np.savetxt(
                os.path.join(debug_dir, "prior_var.txt"),
                [np.nan] if (prior_var is None) else [prior_var],
            )
            np.savetxt(
                os.path.join(debug_dir, "prior_weights.txt"),
                prior_weights if use_prior_causal_prob else [np.nan],
            )

            # create a zipped debug file
            import zipfile

            debug_files = glob.glob(os.path.join(debug_dir, "*.txt"))
            zip_file = os.path.join(debug_dir, "debug.zip")
            zf = zipfile.ZipFile(zip_file, mode="w")
            for debug_file in debug_files:
                zf.write(
                    debug_file,
                    os.path.basename(debug_file),
                    compress_type=zipfile.ZIP_DEFLATED,
                )

        assert self.df_ld.notnull().all().all()
        if residual_var is not None:
            residual_var_init = residual_var

        if hasattr(self.susieR, "susie_suff_stat"):
            logging.info("Using susieR::susie_suff_stat()")
            # susie_suff_stat => susie_rss; this is works for new version of susieR and no diff
            susie_obj = self.susieR.susie_rss(
                bhat=bhat.reshape((m, 1)),
                shat=np.ones((m, 1)),
                R=self.df_ld.values,
                n=self.n,
                L=num_causal_snps,
                scaled_prior_variance=(0.0001 if (prior_var is None) else prior_var),
                estimate_prior_variance=(prior_var is None),
                residual_variance=(
                    self.R_null if (residual_var_init is None) else residual_var_init
                ),
                estimate_residual_variance=(residual_var is None),
                max_iter=susie_max_iter,
                verbose=verbose,
                prior_weights=(
                    prior_weights.reshape((m, 1))
                    if use_prior_causal_prob
                    else self.R_null
                ),
            )
        elif hasattr(self.susieR, "susie_bhat"):
            logging.info("Using susieR::susie_bhat()")
            susie_obj = self.susieR.susie_bhat(
                bhat=bhat.reshape((m, 1)),
                shat=np.ones((m, 1)),
                R=self.df_ld.values,
                n=self.n,
                L=num_causal_snps,
                scaled_prior_variance=(0.0001 if (prior_var is None) else prior_var),
                estimate_prior_variance=(prior_var is None),
                residual_variance=(
                    self.R_null if (residual_var is None) else residual_var
                ),
                estimate_residual_variance=(residual_var is None),
                max_iter=susie_max_iter,
                verbose=verbose,
                prior_weights=(
                    prior_weights.reshape((m, 1))
                    if use_prior_causal_prob
                    else self.R_null
                ),
            )
        else:
            raise NotImplementedError(
                "Only susie_suff_stat() and susie_bhat() are supported. Check your version of susieR"
            )
        susie_time = time.time() - t0
        logging.info("Done in %0.2f seconds" % (susie_time))

        # extract pip and beta_mean
        pip = np.array(self.susieR.susie_get_pip(susie_obj))
        beta_mean = np.array(self.susieR.coef_susie(susie_obj)[1:])
        assert np.allclose(
            beta_mean,
            np.sum(
                np.array(susie_obj.rx2("mu")) * np.array(susie_obj.rx2("alpha")), axis=0
            )
            / np.array(susie_obj.rx2("X_column_scale_factors")),
        )

        # compute the posterior mean of beta^2
        s_alpha = np.array(susie_obj.rx2("alpha"))
        s_mu = np.array(susie_obj.rx2("mu"))
        s_mu2 = np.array(susie_obj.rx2("mu2"))
        s_X_column_scale_factors = np.array(susie_obj.rx2("X_column_scale_factors"))
        beta_var = np.sum(s_alpha * s_mu2 - (s_alpha * s_mu) ** 2, axis=0) / (
            s_X_column_scale_factors**2
        )
        assert np.all(beta_var >= 0)

        # create output df
        df_susie = self.df_sumstats_locus.copy()
        df_susie["PIP"] = pip
        df_susie["BETA_MEAN"] = beta_mean
        # flip back the finemap BETA, as the alleles are in original order
        df_susie.loc[is_flipped, "BETA_MEAN"] *= -1
        df_susie["BETA_SD"] = np.sqrt(beta_var)

        # add distance from center
        start = df_susie["BP"].min()
        end = df_susie["BP"].max()
        middle = (start + end) // 2
        df_susie["DISTANCE_FROM_CENTER"] = np.abs(df_susie["BP"] - middle)

        # mark causal sets
        import rpy2

        logging.info("Using rpy2 version %s" % (rpy2.__version__))
        if Version(rpy2.__version__) >= Version("3.5.9"):
            snames = (susie_obj.names).tolist()
            self.susie_dict = {
                key: np.array(susie_obj.rx2(key), dtype=object) for key in snames
            }
        else:
            self.susie_dict = {
                key: np.array(susie_obj.rx2(key), dtype=object)
                for key in list(susie_obj.names)
            }
        df_susie["CREDIBLE_SET"] = 0
        susie_sets = self.susie_dict["sets"][0]
        # if type(susie_sets) != self.RNULLType:
        try:
            for set_i, susie_set in enumerate(susie_sets):
                is_in_set = np.zeros(df_susie.shape[0], dtype=bool)
                is_in_set[np.array(susie_set) - 1] = True
                is_in_set[df_susie["CREDIBLE_SET"] > 0] = False
                df_susie.loc[is_in_set, "CREDIBLE_SET"] = set_i + 1
        except TypeError:
            pass

        # save SuSiE object if requested
        if susie_outfile is not None:
            from rpy2.robjects.packages import importr

            R_base = importr(
                "base",
                robject_translations={
                    "print.me": "print_dot_me",
                    "print_me": "print_uscore_me",
                },
            )
            R_base.saveRDS(susie_obj, file=susie_outfile)
            logging.info("Saved SuSiE object to RDS file: %s" % (susie_outfile))

        # delete the LD file if needed
        if delete_ld_files_on_exit:
            ld_file_dir = os.path.dirname(ld_file)
            if os.path.exists(ld_file_dir):
                shutil.rmtree(ld_file_dir)

        return df_susie
