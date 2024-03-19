from finemap_tools.reader.gwas.sumstats import load_sumstats
import numpy as np
import pandas as pd
import os
import time
import scipy.stats as stats
import logging

from tqdm import tqdm
import tempfile
import glob
from finemap_tools.others.polyfun.polyfun_utils import set_snpid_index
from finemap_tools.others.polyfun.ldstore.bcor import bcor
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer
from finemap_tools.utils import run_executable
from scipy import sparse
from finemap_tools.reader.ld import read_ld_from_file
# from polyfun import configure_logger, check_package_versions

from pathlib import Path


def save_ld_to_npz(ld_arr, df_ld_snps, npz_file):

    assert npz_file.endswith(".npz")
    logging.info("Saving LD file %s" % (npz_file))
    t0 = time.time()

    # save meta file
    meta_file = npz_file[:-4] + ".gz"
    df_ld_snps.to_csv(meta_file, sep="\t", index=False)

    # save .npz file
    R = np.tril(ld_arr).astype(np.float64)
    np.fill_diagonal(R, np.diag(R) / 2.0)
    R = sparse.coo_matrix(R)
    sparse.save_npz(npz_file, R, compressed=True)
    logging.info("Done in %0.2f seconds" % (time.time() - t0))

class Fine_Mapping(object):
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

        # check that data is valid
        if genotypes_file is not None:
            if genotypes_file.endswith(".bgen"):
                if sample_file is None:
                    raise IOError("sample-file must be provided with a bgen file")

        df_sumstats = load_sumstats(
            sumstats_file,
            chr_num,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
        )

        # save class members
        self.genotypes_file = genotypes_file
        self.n = n
        self.sample_file = sample_file
        self.df_sumstats = df_sumstats
        self.incl_samples = incl_samples
        self.ldstore_exe = ldstore_exe
        self.cache_dir = cache_dir
        self.cache_format = cache_format
        self.n_threads = n_threads
        self.chr = chr_num
        self.memory = memory
        self.allow_swapped_indel_alleles = allow_swapped_indel_alleles

    def sync_ld_sumstats(self, ld_arr, df_ld_snps, allow_missing=False):
        df_ld_snps = set_snpid_index(
            df_ld_snps, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )

        if ld_arr is None:
            df_ld = pd.DataFrame(
                np.zeros(len(df_ld_snps.index), dtype=np.int64),
                index=df_ld_snps.index,
                columns=["dummy"],
            )
        else:
            assert ld_arr.shape[0] == df_ld_snps.shape[0]
            assert ld_arr.shape[0] == ld_arr.shape[1]
            df_ld = pd.DataFrame(
                ld_arr, index=df_ld_snps.index, columns=df_ld_snps.index
            )
        # TODO: rm some LD with NaNs

        have_na_snp = df_ld[df_ld.isna().sum() >= 1].index.tolist()
        logging.info(
            f"remove {len(have_na_snp)} SNPs with NaNs in LD, this should be think carefully."
        )
        df_ld = df_ld.drop(index=have_na_snp, columns=have_na_snp)

        # make sure that all SNPs in the sumstats file are in the LD file
        snps_in_ld_file = self.df_sumstats_locus.index.isin(df_ld.index)
        if not np.all(snps_in_ld_file):
            # Could the missing SNPs be due to mismatched indel alleles?
            df_sumstats_missing = self.df_sumstats_locus.loc[~snps_in_ld_file]
            num_missing_is_indel = np.sum(
                (df_sumstats_missing["A1"].str.len() > 1)
                | (df_sumstats_missing["A2"].str.len() > 1)
            )
            if allow_missing:
                num_missing = np.sum(~snps_in_ld_file)
                logging.warning(
                    "%d variants with sumstats were not found in the LD file and will be omitted (please note that this may lead to false positives if the omitted SNPs are causal!)"
                    % (num_missing)
                )
                if num_missing_is_indel > 0 and not self.allow_swapped_indel_alleles:
                    logging.warning(
                        "%d of the missing variants were indels. Check that the allele order (A1/A2) matches between the sumstats and the LD file. Also see the flag --allow-swapped-indel-alleles"
                        % (num_missing_is_indel)
                    )
                self.df_sumstats_locus = self.df_sumstats_locus.loc[snps_in_ld_file]
                assert np.all(self.df_sumstats_locus.index.isin(df_ld.index))
            else:
                error_msg = (
                    "not all SNPs in the sumstats file were found in the LD matrix!"
                    " You could drop the missing SNPs with the flag --allow-missing, but please note that"
                    " these omitted SNPs may be causal, in which case you may get false positive results..."
                    " If there should be no missing SNPs (e.g. you are using insample LD), see the flag --allow-swapped-indel-alleles"
                )
                raise ValueError(error_msg)

        # make sure that df_sumstats_locus is not empty
        assert (
            self.df_sumstats_locus.shape[0] > 0
        ), "no SNPs found in df_sumstats_locus. Please double-check that the SNPs in the LD file and in the sumstats file have the exact same positions"

        # filter LD to only SNPs found in the sumstats file
        assert not np.any(self.df_sumstats_locus.index.duplicated())
        if df_ld.shape[0] != self.df_sumstats_locus.shape[0] or np.any(
            df_ld.index != self.df_sumstats_locus.index
        ):
            if ld_arr is None:
                df_ld = df_ld.loc[self.df_sumstats_locus.index]
            else:
                df_ld = df_ld.loc[
                    self.df_sumstats_locus.index, self.df_sumstats_locus.index
                ]
            df_ld_snps = df_ld_snps.loc[df_ld.index]

        # do a final verification that we're synced
        assert np.all(df_ld.index == self.df_sumstats_locus.index)
        assert np.all(df_ld_snps.index == self.df_sumstats_locus.index)

        # add leading zero to sumstats CHR column if needed
        if np.any(df_ld_snps["CHR"].astype(str).str.startswith("0")):
            self.df_sumstats_locus = self.df_sumstats_locus.copy()
            self.df_sumstats_locus["CHR"] = self.df_sumstats_locus["CHR"].astype(str)
            is_1digit = self.df_sumstats_locus["CHR"].str.len() == 1
            self.df_sumstats_locus.loc[is_1digit, "CHR"] = (
                "0" + self.df_sumstats_locus.loc[is_1digit, "CHR"]
            )

        # update self.df_ld
        self.df_ld = df_ld
        self.df_ld_snps = df_ld_snps
        assert self.df_ld.notnull().all().all()

    def find_cached_ld_file(self, locus_start, locus_end, need_bcor=False):

        # if there's no cache dir, return None
        if self.cache_dir is None:
            return None

        if self.incl_samples is None:
            fname_pattern = "%s.%d" % (os.path.basename(self.genotypes_file), self.chr)
        else:
            fname_pattern = "%s.%s.%d" % (
                os.path.basename(self.genotypes_file),
                os.path.basename(self.incl_samples),
                self.chr,
            )

        # search for suitable LD files
        bcor_files = glob.glob(os.path.join(self.cache_dir, fname_pattern + "*.bcor"))
        npz_files = glob.glob(os.path.join(self.cache_dir, fname_pattern + "*.npz"))
        if need_bcor:
            ld_files = bcor_files + npz_files
        else:
            ld_files = npz_files + bcor_files

        for ld_file in ld_files:
            if os.stat(ld_file).st_size == 0:
                os.remove(ld_file)
                continue
            ld_basename = os.path.basename(ld_file)
            bp1 = int(ld_basename.split(".")[-3])
            bp2 = int(ld_basename.split(".")[-2])
            assert bp1 < bp2
            if (bp1 > locus_start) or (bp2 < locus_end):
                continue

            # get the list of SNPs in the LD file
            if ld_file.endswith(".npz"):
                meta_file = ld_file[:-4] + ".gz"
                if not os.path.exists(meta_file):
                    continue
                df_ld_snps = pd.read_table(
                    meta_file,
                    sep="\s+",
                    usecols=["A1", "A2", "BP", "CHR", "SNP"],
                )
            elif ld_file.endswith(".bcor"):
                bcor_obj = bcor(ld_file)
                df_ld_snps = bcor_obj.getMeta()
                del bcor_obj
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
            else:
                raise IOError("unknown file extension")

            df_ld_snps = set_snpid_index(
                df_ld_snps, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
            )

            # make sure that the LD file includes data for all the SNPs in the locus
            if not np.all(self.df_sumstats_locus.index.isin(df_ld_snps.index)):
                logging.warning(
                    "The available cached LD file was ignored because it does not contain data for all the SNPs in the locus"
                )
                continue

            # if we got here than we found a suitable d file
            logging.info(
                "Found a cached LD file containing all SNPs with sumstats in chromosome %d BP %d-%d: %s"
                % (self.chr, locus_start, locus_end, ld_file)
            )
            return ld_file

    def get_ld_output_file_prefix(self, locus_start, locus_end, output_dir=None):
        if self.cache_dir is None:
            if output_dir is None:
                output_dir = tempfile.mkdtemp()
            output_prefix = os.path.join(output_dir, "ld")
        else:
            if self.incl_samples is None:
                output_prefix = os.path.join(
                    self.cache_dir,
                    "%s.%d.%d.%d"
                    % (
                        os.path.basename(self.genotypes_file),
                        self.chr,
                        locus_start,
                        locus_end,
                    ),
                )
            else:
                output_prefix = os.path.join(
                    self.cache_dir,
                    "%s.%s.%d.%d.%d"
                    % (
                        os.path.basename(self.genotypes_file),
                        os.path.basename(self.incl_samples),
                        self.chr,
                        locus_start,
                        locus_end,
                    ),
                )

        return output_prefix

    def compute_ld_bgen(self, locus_start, locus_end, verbose=False):

        # create df_z
        df_z = self.df_sumstats_locus[["SNP", "CHR", "BP", "A1", "A2"]].copy()

        try:
            import bgen
        except (ImportError, ModuleNotFoundError):
            raise ValueError(
                '\n\nPlease install the bgen package (using "pip install bgen")'
            )
        from bgen.reader import BgenFile

        bfile = BgenFile(self.genotypes_file)
        bgen_chromosomes = bfile.chroms()
        if bgen_chromosomes[0].startswith("0"):
            df_z["CHR"] = "0" + df_z["CHR"].astype(str)

        # sync the order of the alleles between the sumstats and the bgen file
        list_bgen = []
        # rsids = bfile.rsids()
        # small change reduces the time for bgen processing
        # the previous implementation would iterate through all the SNPs in the bgen file
        # this implementation loops over just the snps in the locus
        rsids = bfile.fetch(self.chr, locus_start, locus_end)
        for snp_i, rsid in enumerate(rsids):
            #             if rsid not in df_z['SNP'].values: continue
            #             snp_alleles = bfile[snp_i].alleles
            #             snp_chrom = bfile[snp_i].chrom
            #             snp_pos = bfile[snp_i].pos
            if rsid.rsid not in df_z["SNP"].values:
                continue  # NOTE: so this indeed take the intersection of the two sets
            snp_alleles = rsid.alleles
            snp_chrom = rsid.chrom
            snp_pos = rsid.pos
            rsid = rsid.rsid  # NOTE: this is the change
            assert (
                len(snp_alleles) == 2
            ), "cannot handle SNPs with more than two alleles"
            df_snp = df_z.query('SNP == "%s"' % (rsid))

            assert df_snp.shape[0] == 1
            a1, a2 = df_snp["A1"].iloc[0], df_snp["A2"].iloc[0]
            snp_a1, snp_a2 = snp_alleles[0], snp_alleles[1]
            if set([a1, a2]) != set([snp_a1, snp_a2]):
                raise ValueError(
                    "The alleles for SNP %s are different in the sumstats and in the bgen file:\n \
                                 bgen:     A1=%s  A2=%s\n \
                                 sumstats: A1=%s  A2=%s \
                                "
                    % (rsid, snp_alleles[0], snp_alleles[1], a1, a2)
                )
            d = {
                "SNP": rsid,
                "CHR": snp_chrom,
                "BP": snp_pos,
                "A1": snp_a1,
                "A2": snp_a2,
            }
            list_bgen.append(d)
        df_bgen = pd.DataFrame(list_bgen)
        df_bgen = set_snpid_index(
            df_bgen, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )
        df_z = set_snpid_index(
            df_z, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )
        df_z = df_z[[]].merge(df_bgen, left_index=True, right_index=True)
        df_z = df_z[["SNP", "CHR", "BP", "A1", "A2"]]

        # rename columns
        df_z.rename(
            columns={
                "SNP": "rsid",
                "CHR": "chromosome",
                "BP": "position",
                "A1": "allele1",
                "A2": "allele2",
            },
            inplace=True,
        )

        # Create LDstore input files
        temp_dir = tempfile.mkdtemp()
        incl_file = os.path.join(
            temp_dir, "incl.incl"
        )  # NOTE: remove this line if not needed
        master_file = os.path.join(temp_dir, "master.master")
        z_file = os.path.join(
            temp_dir, "chr%s.%s_%s.z" % (self.chr, locus_start, locus_end)
        )
        dose_file = os.path.join(temp_dir, "dosages.bdose")
        df_z.to_csv(z_file, sep=" ", index=False)

        # find number of samples
        if (
            self.incl_samples is None
        ):  # NOTE: looks prety  bad as this is only the reason of wheather there is a header in it
            num_samples = pd.read_table(self.sample_file).shape[0] - 1
        else:
            num_samples = pd.read_table(self.incl_samples, header=None).shape[0]

        # get output file name
        bcor_file = os.path.join(
            temp_dir, "chr%s.%s_%s.bcor" % (self.chr, locus_start, locus_end)
        )
        ## parse for cache_format
        if self.cache_format == "npz":  #
            ld_file = os.path.join(
                temp_dir, "chr%s.%s_%s.ld" % (self.chr, locus_start, locus_end)
            )

        elif self.cache_format == "bcor":
            bcor_file = (
                self.get_ld_output_file_prefix(locus_start, locus_end, temp_dir)
                + ".bcor"
            )
        else:
            raise ValueError(f"unknown cache format {self.cache_format}")

        # Create LDstore master file
        df_master = pd.DataFrame(
            columns=["z", "bgen", "bgi", "bcor", "dose", "sample", "n_samples"]
        )
        df_master["z"] = [z_file]
        df_master["bgen"] = [self.genotypes_file]
        df_master["bgi"] = [self.genotypes_file + ".bgi"]
        df_master["bcor"] = [bcor_file]
        df_master["ld"] = [ld_file]
        df_master["bdose"] = [dose_file]
        df_master["sample"] = [self.sample_file]
        df_master["n_samples"] = num_samples
        if self.incl_samples is not None:
            df_master["incl"] = self.incl_samples
        df_master.to_csv(master_file, sep=";", header=True, index=False)

        # run LDstore
        ldstore_cmd = [
            self.ldstore_exe,
            "--in-files",
            master_file,
            # "--write-bcor",
            # "--write-text",
            "--write-bdose",
            "--bdose-version",
            "1.0",
        ]  # TODO: maybe for checking big files or for bdose 1.1

        if self.cache_format == "npz":
            ldstore_cmd += ["--write-text"]
        elif self.cache_format == "bcor":
            ldstore_cmd += ["--write-bcor"]

        if self.memory is not None:
            ldstore_cmd += ["--memory", str(self.memory)]
        if self.n_threads is not None:
            ldstore_cmd += ["--n-threads", str(self.n_threads)]
        run_executable(
            ldstore_cmd,
            "LDStore",
            measure_time=True,
            show_output=verbose,
            show_command=verbose,
        )

        if self.cache_format == "bcor":
            if not os.path.exists(bcor_file):
                raise IOError("Could not find output BCOR file")
            return bcor_file
        elif (
            self.cache_format == "npz"
        ):  # load txt file by np and return ld_arr and df_z
            if not os.path.exists(ld_file):
                raise IOError("Could not find output LD file")

            ld_arr = np.loadtxt(ld_file)
            assert ld_arr.shape[0] == ld_arr.shape[1] == df_z.shape[0]
            df_ld_snps = df_z.rename(
                columns={
                    "allele1": "A1",
                    "allele2": "A2",
                    "position": "BP",
                    "chromosome": "CHR",
                    "rsid": "SNP",
                },
                inplace=False,
                errors="ignore",
            )
            # save_ld_to_npz(ld_arr, df_z, npz_file) # NOTE: after this function will auto save to npz file, so only need to return ld_arr, df_z
            # return
            return ld_arr, df_ld_snps

    def read_plink_genotypes(self, bed):
        X = bed.compute().astype(np.float64)
        if np.any(np.isnan(X)): # is this ok ????
            imp = SimpleImputer(missing_values=np.nan, strategy="mean", copy=False)
            imp.fit(X)
            X = imp.transform(X)
        X -= X.mean(axis=0)
        assert not np.any(np.isnan(X))
        is_polymorphic = X.std(axis=0) > 0
        X[:, is_polymorphic] /= X[:, is_polymorphic].std(axis=0)
        return X

    def compute_ld_plink_pgen(self, locus_start, locus_end, verbose):
        from finemap_tools.plink import plink2_cal_LD
        from finemap_tools.reader.plink import read_pvar

        logging.info(
            f"Computing LD from pgen fileset {self.genotypes_file} chromosome {self.chr} region {locus_start}-{locus_end}"
        )
        t0 = time.time()

        ld_snp_df = read_pvar(Path(self.genotypes_file).with_suffix(".pvar")).rename(
            columns={
                "#CHROM": "CHR",
                "ID": "SNP",
                "POS": "BP",
                "REF": "A2",  # REF as A2 is to make sure the A1 is the minor allele and consistent with the sumstats later
                "ALT": "A1",
            }
        )
        ld_snp_df = set_snpid_index(
            ld_snp_df, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )
        logging.info(f"read {ld_snp_df.shape[0]} SNPs from pvar file")

        # used sumstats to include the snp in analysis
        df_z = self.df_sumstats_locus[["SNP", "CHR", "BP", "A1", "A2"]].copy()
        df_z = set_snpid_index(
            df_z, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )
        df_z = df_z[[]].merge(ld_snp_df, left_index=True, right_index=True)
        ld_snp_df = df_z[["SNP", "CHR", "BP", "A1", "A2"]]
        del df_z

        logging.info(
            f"only keep {ld_snp_df.shape[0]} SNPs in the sumstats file to cal LD to avoid the long time of cal LD"
        )

        # get_ld
        with tempfile.TemporaryDirectory() as tmpdirname:
            print("created temporary directory", tmpdirname)
            ld_df = plink2_cal_LD(
                pgen=self.genotypes_file,
                snplist=ld_snp_df["SNP"].tolist(),
                outSuffix=tmpdirname + "/test",
                thread=self.n_threads,
                memory=self.memory * 1024,
            )

        assert ld_df.shape[0] == ld_snp_df.shape[0]
        assert (
            ld_df.index.tolist() == ld_snp_df["SNP"].tolist() == ld_df.columns.tolist()
        )

        return ld_df.values, ld_snp_df

    def compute_ld_plink(self, locus_start, locus_end, verbose):
        logging.info(
            "Computing LD from plink fileset %s chromosome %s region %s-%s"
            % (self.genotypes_file, self.chr, locus_start, locus_end)
        )
        t0 = time.time()

        # read the plink file
        df_bim, df_fam, bed = read_plink(self.genotypes_file)
        df_bim.rename(
            columns={"snp": "SNP", "pos": "BP", "chrom": "CHR", "a0": "A2", "a1": "A1"},
            inplace=True,
        )
        df_bim["A1"] = df_bim["A1"].astype("str")
        df_bim["A2"] = df_bim["A2"].astype("str")
        df_bim["CHR"] = df_bim["CHR"].astype(np.int64)
        del df_bim["i"]
        del df_bim["cm"]
        bed = bed.T

        # zoom in on target locus
        is_snp_in_region = (df_bim["BP"].between(locus_start, locus_end)) & (
            df_bim["CHR"] == self.chr
        )
        df_bim = df_bim.loc[is_snp_in_region]
        df_ld_snps = df_bim
        bed = bed[:, is_snp_in_region.values]
        assert bed.shape[1] > 0, "No SNPs found in the target region"

        # compute chunk size, using the formula MEM = bed.shape[0] * chunk_size * 4 / 2**30
        if self.memory is None:
            mem_limit = 1
        else:
            mem_limit = self.memory
        chunk_size = np.int64(
            (np.float64(mem_limit) * 0.8) / bed.shape[0] / 4 * (2**30)
        )
        if chunk_size == 0:
            chunk_size = 1
        if chunk_size > bed.shape[1]:
            chunk_size = bed.shape[1]
        num_chunks = np.int64(np.ceil(bed.shape[1] / chunk_size))
        if num_chunks > 1:
            assert chunk_size * (num_chunks - 2) < bed.shape[1] - 1
        if chunk_size * (num_chunks - 1) >= bed.shape[1]:
            num_chunks -= 1

        # compute LD in chunks
        logging.info(
            "Found %d SNPs in target region. Computing LD in %d chunks..."
            % (bed.shape[1], num_chunks)
        )
        ld_arr = np.empty((bed.shape[1], bed.shape[1]), dtype=np.float64)
        for chunk_i in tqdm(range(num_chunks)):
            chunk_i_start = chunk_i * chunk_size
            chunk_i_end = np.minimum(chunk_i_start + chunk_size, bed.shape[1])
            X_i = self.read_plink_genotypes(bed[:, chunk_i_start:chunk_i_end])
            ld_arr[chunk_i_start:chunk_i_end, chunk_i_start:chunk_i_end] = (
                X_i.T.dot(X_i) / X_i.shape[0]
            )
            for chunk_j in range(chunk_i + 1, num_chunks):
                chunk_j_start = chunk_j * chunk_size
                chunk_j_end = np.minimum(chunk_j_start + chunk_size, bed.shape[1])
                X_j = self.read_plink_genotypes(bed[:, chunk_j_start:chunk_j_end])
                ld_arr[chunk_i_start:chunk_i_end, chunk_j_start:chunk_j_end] = (
                    X_i.T.dot(X_j) / X_i.shape[0]
                )
                ld_arr[chunk_j_start:chunk_j_end, chunk_i_start:chunk_i_end] = ld_arr[
                    chunk_i_start:chunk_i_end, chunk_j_start:chunk_j_end
                ].T
        ld_arr = np.nan_to_num(ld_arr, copy=False)
        ld_diag = np.diag(ld_arr).copy()
        if np.any(np.isclose(ld_diag, 0.0)):
            ld_diag[np.isclose(ld_diag, 0.0)] = 1.0
            np.fill_diagonal(ld_arr, ld_diag)

        logging.info("Done in %0.2f seconds" % (time.time() - t0))
        return ld_arr, df_ld_snps

    def set_locus(self, locus_start, locus_end):

        # update self.df_sumstats_locus
        self.df_sumstats_locus = self.df_sumstats.query(
            "%d <= BP <= %d" % (locus_start, locus_end)
        )
        num_snps = self.df_sumstats_locus.shape[0]
        if num_snps < 2:
            raise ValueError(
                "%d SNP(s) found in sumstats file in the BP range %d-%d"
                % (num_snps, locus_start, locus_end)
            )

    def get_ld_data(self, locus_start, locus_end, need_bcor=False, verbose=False):

        ld_arr, df_ld_snps, ld_file = None, None, None

        # check if we already have a suitable LD file in the cache dir
        ld_file = self.find_cached_ld_file(locus_start, locus_end, need_bcor=need_bcor)

        # compute LD if we couldn't find a suitable LD file
        if ld_file is None:
            if self.genotypes_file.endswith(".bgen"):  # this won't return None
                if not os.path.exists(self.genotypes_file):
                    raise IOError("%s doesn't exist" % (self.genotypes_file))
                if self.cache_format == "bcor":
                    ld_file = self.compute_ld_bgen(
                        locus_start, locus_end, verbose=verbose
                    )
                elif self.cache_format == "npz":
                    ld_arr, df_ld_snps = self.compute_ld_bgen(
                        locus_start, locus_end, verbose=verbose
                    )
                # ld_file = self.compute_ld_bgen(locus_start, locus_end, verbose=verbose)
            elif os.path.exists(self.genotypes_file + ".bed"):
                ld_arr, df_ld_snps = self.compute_ld_plink(
                    locus_start, locus_end, verbose=verbose
                )
            elif os.path.exists(self.genotypes_file + ".pgen"):
                ld_arr, df_ld_snps = self.compute_ld_plink_pgen(
                    locus_start, locus_end, verbose=verbose
                )
            else:
                raise ValueError(
                    "no suitable file found for: %s" % (self.genotypes_file)
                )

        # arrange the LD data
        assert ld_file is None or (ld_arr is None and df_ld_snps is None)

        # if there is no LD file, return the LD data directly

        if ld_file is None:  #  NOTE: this is not possiblely be None as the code
            # cache output if possible
            if self.cache_dir is not None:
                npz_file = (
                    self.get_ld_output_file_prefix(locus_start, locus_end) + ".npz"
                )
                save_ld_to_npz(ld_arr, df_ld_snps, npz_file)
            return ld_arr, df_ld_snps
        # if we have an LD file, return it if it's a bcor and we want a bcor, or return the LD data directly otherwise
        else:  # NOTE: this code is quiet strange
            if ld_file.endswith(".bcor"):

                if need_bcor:
                    return ld_file
                else:
                    ld_arr, df_ld_snps = read_ld_from_file(ld_file)

                    # cache output if possible
                    if self.cache_dir is not None and ld_file.endswith(".bcor"):
                        npz_file = (
                            self.get_ld_output_file_prefix(locus_start, locus_end)
                            + ".npz"
                        )
                        save_ld_to_npz(ld_arr, df_ld_snps, npz_file)
                    return ld_arr, df_ld_snps

            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                return ld_arr, df_ld_snps

    def finemap(self):
        raise NotImplementedError()

    def estimate_h2_hess(self, prop_keep=0.005, R_cutoff=0.99, pvalue_bound=None):
        """
        prop_keep:  Proportion of SNPs to use in the estimation (only the ones with the smallest p-values)
        R_cutoff: Exclude one of each pair of SNPs with with magnitude of correlation greater than this value
        pvalue_bound: An upper bound on the p-value cutoff (i.e., SNPs with P greater than this cutoff will never be used in the estimation)

        The modified HESS equation implemented below is

        $$ \frac{n \alpha R^{-1} \alpha - m}{n} = \alpha R^{-1} \alpha - \frac{m}{n} $$

        where $\alpha$ is a vector of marginal effect size estimates for $m$ standardized SNPs,
        $R$ is a matrix of summary LD information, and $n$ is the sample size.

        This is a biased estimator (denominator of $n$) with a smaller estimation variance compared
        with the unbiased estimator (denominator of $n - m$) used in the original HESS publication
        (Shi et al., 2014; https://doi.org/10.1016/j.ajhg.2016.05.013).
        """

        # keep only potential causal SNPs
        pvalue_cutoff = self.df_sumstats_locus["P"].quantile(prop_keep)
        if pvalue_cutoff == 0:
            pvalue_cutoff = np.min(self.df_sumstats_locus["P"].loc[lambda p: p > 0])
        if pvalue_bound is not None and pvalue_cutoff > pvalue_bound:
            pvalue_cutoff = pvalue_bound
        is_potential_csnp = self.df_sumstats_locus["P"].values < pvalue_cutoff
        if np.any(is_potential_csnp):
            R_pot_csnp = self.df_ld.loc[is_potential_csnp, is_potential_csnp].values
        else:
            return 0

        # take a maximally independent subset
        np.fill_diagonal(R_pot_csnp, 0)
        import networkx as nx

        G = nx.from_numpy_array(np.abs(R_pot_csnp) > R_cutoff)
        np.fill_diagonal(R_pot_csnp, 1)
        inds = np.sort(nx.maximal_independent_set(G))

        # estimate h2 using HESS
        R_subset = R_pot_csnp[np.ix_(inds, inds)]
        alpha_subset = self.df_sumstats_locus.loc[is_potential_csnp, "Z"].iloc[
            inds
        ].values / np.sqrt(self.n)
        h2_hess = (
            alpha_subset.dot(np.linalg.solve(R_subset, alpha_subset))
            - R_subset.shape[0] / self.n
        )

        return h2_hess

    def estimate_h2_hess_wrapper(
        self, prop_keep=0.005, R_cutoff=0.99, min_h2=None, num_samples=100
    ):
        """
        prop_keep:  Proprtion of SNPs to use in the estimation (only the ones with the smallest p-values)
        R_cutoff: Exclude one of each pair of SNPs with with magnitude of correlation greater than this value
        min_h2: Exclude SNPs that tag less than this amount of heritability
        num_samples: Number of random samples of indepdendent SNPs to draw
        """

        if min_h2 is None:
            pvalue_bound = None
        else:
            assert (
                min_h2 > 0 and min_h2 < 1
            ), "The minimum proportion of heritability to exclude SNPs from HESS estimation must be between 0 and 1"
            pvalue_bound = stats.chi2(1).sf(min_h2 * self.n)

        assert num_samples > 0, "Number of random samples must be a positive integer"

        h2_hess_list = [
            self.estimate_h2_hess(
                prop_keep=prop_keep, R_cutoff=R_cutoff, pvalue_bound=pvalue_bound
            )
            for try_num in range(num_samples)
        ]
        h2_hess = np.mean(h2_hess_list)
        return h2_hess
