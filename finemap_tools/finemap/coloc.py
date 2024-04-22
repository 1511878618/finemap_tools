# from .base import Fine_Mapping
from finemap_tools.reader.ld import read_ld_from_file
from finemap_tools.others.polyfun.polyfun_utils import set_snpid_index
import numpy as np
from finemap_tools.reader.gwas.sumstats import load_sumstats
from finemap_tools.finemap.base import Fine_Mapping
import logging
from finemap_tools.utils import add_ID
from pathlib import Path


default_coloc_map = {
    "snp": "snp",
    "chr": "chrom",
    "position": "pos",
    "beta": "beta",
    "varbeta": None,
}


necessary_check_cols = ["beta", "varbeta"]  # drop na by this col


## load data
def sumstats2coloc(sumstats, map_dict=None):
    # sumstats = sumstats.copy()
    polyfun_cols = list(default_coloc_map.keys())
    used_polyfun_map_dict = default_coloc_map.copy()
    if map_dict:
        used_polyfun_map_dict.update(map_dict)

    for col in polyfun_cols:
        current_check_col = used_polyfun_map_dict.get(col, None)

        if (current_check_col not in sumstats.columns) | (current_check_col is None):
            if col == "varbeta":
                # assert beta and se in map_dict
                assert (
                    "sebeta" in used_polyfun_map_dict.keys()
                ), "sebeta not in map_dict"
                sumstats["varbeta"] = sumstats[used_polyfun_map_dict["sebeta"]] ** 2

                used_polyfun_map_dict["varbeta"] = "varbeta"
            elif col == "snp":
                sumstats["snp"] = add_ID(
                    sumstats,
                    col_list=[
                        used_polyfun_map_dict["chr"],
                        used_polyfun_map_dict["position"],
                        used_polyfun_map_dict["ref"],
                        used_polyfun_map_dict["alt"],
                    ],
                )
            else:

                raise ValueError(f"{col} and value {current_check_col} not in sumstats")

    rename_dict = {v: k for k, v in used_polyfun_map_dict.items() if v is not None}
    sumstats.rename(columns=rename_dict, inplace=True)
    return sumstats[
        list(rename_dict.values())
        + [col for col in sumstats.columns if col not in rename_dict.values()]
    ]


def load_summstats_coloc(
    filepath,
    map_dict=None,
    chr_num=None,
    start=None,
    end=None,
    allow_swapped_indel_alleles=False,
    sumstats_format=None,
):

    if filepath.endswith(".gz") and Path(filepath + ".tbi").exists():
        if sumstats_format == "GTExV8":
            logging.info(f"loading GTExV8 format sumstats file with {filepath}......")
            from finemap_tools.reader.eQTL.GTEx import GTEx_tabix_reader

            df_sumstats = GTEx_tabix_reader(
                tabix_file=filepath, region=f"{chr_num}:{start}-{end}"
            )
        else:
            logging.info(
                f"loading with normal tabix reader to loading sumstats file with {filepath}......"
            )
            from finemap_tools.reader import tabix_reader

            df_sumstats = tabix_reader(filepath, region=f"{chr_num}:{start}-{end}")
    else:
        import pandas as pd
        import numpy as np

        try:
            df_sumstats = pd.read_table(filepath, sep="\s+")
        except:
            raise ValueError(f"cannot read file {filepath} with pd.read_table")

        logging.info(
            f"loading sumstats file with {filepath}, use normal pd.read_table with rows loaded: {df_sumstats.shape[0]}; this should be only used with user split sumstats file and not recommended for large sumstats file"
        )

    logging.info(
        f"the header of loaded file is {df_sumstats.columns}, if it is not the correct header, please check the file or with correct format specified."
    )

    if df_sumstats.shape[0] == 0:
        raise IOError(
            f"sumstats file does not include any SNPs in chromosome {chr_num} from {start} to {end}"
        )

    df_sumstats = sumstats2coloc(df_sumstats, map_dict=map_dict)

    from finemap_tools.snpfilter import filter_pipline

    logging.info(
        f"filtering SNP by finemap_tools with {df_sumstats.shape[0]} SNP at begining........"
    )
    df_sumstats = filter_pipline(sumstats=df_sumstats, id_col="snp")

    logging.info(f"after filtering, left {df_sumstats.shape[0]} SNP")

    is_na_nums_sumstats = df_sumstats[necessary_check_cols].isna().any().sum()
    logging.info(f"na numbers in sumstats {Path(filepath).name}: {is_na_nums_sumstats}")

    df_sumstats = df_sumstats.dropna(subset=necessary_check_cols)
    if df_sumstats.shape[0] == 0:
        raise ValueError(f"no valid rows in {Path(filepath).name}")

    return df_sumstats


gwas_coloc_map_dict = {
    "pheweb": {
        "snp": None,
        "chr": "chrom",
        "position": "pos",
        "beta": "beta",
        "varbeta": None,
        "sebeta": "sebeta",
        "alt": "alt",
        "ref": "ref",
        # "A1": "alt",
        # "A2": "ref",
    },
    "polyfun": {
        "SNP": "snp",
        "CHR": "chrom",
        "BP": "pos",
        "A1": "alt",
        "A2": "ref",
    },
}

eQTL_coloc_map_dict = {
    "GTExV8": {
        "snp": None,
        "chr": "chromosome",
        "position": "position",
        "beta": "beta",
        "sebeta": "se",
        "alt": "alt",
        "ref": "ref",
        # "A1": "alt",
        # "A2": "ref",
    }
}

## 寻找需要的列，根据map_dict
## 保留其他列


class coloc_Wrapper(Fine_Mapping):
    def __init__(
        self,
        sumstats_file,
        sumstats_file_format,
        sumstats_file2,
        sumstats_file2_format,
        n1,
        n2,
        type1,
        type2,
        sdY1,
        sdY2,
        chr_num,
        start,
        end,
        sumstats_file_suffix="gwas",
        sumstats_file2_suffix="eQTL",
        ldstore_exe=None,
        genotypes_file=None,
        sample_file=None,
        incl_samples=None,
        cache_dir=None,
        cache_format=None,
        n_threads=None,
        memory=None,
        allow_swapped_indel_alleles=False,  # this do not work currently
    ):

        # check that data is valid
        if genotypes_file is not None:
            if genotypes_file.endswith(".bgen"):
                if sample_file is None:
                    raise IOError("sample-file must be provided with a bgen file")

        df_sumstats = load_sumstats(
            sumstats_file=sumstats_file,
            chr_num=chr_num,
            start=start,
            end=end,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
            sumstats_format=sumstats_file_format,
        )

        if sumstats_file2_format != "pheweb":
            df_sumstats2 = load_summstats_coloc(
                sumstats_file2,
                map_dict=eQTL_coloc_map_dict.get(sumstats_file2_format, None),
                chr_num=chr_num,
                start=start,
                end=end,
                allow_swapped_indel_alleles=allow_swapped_indel_alleles,
                sumstats_format=sumstats_file2_format,
            )
        else:
            # this is a normal gwas dataset not any QTLs with many genes in one file
            df_sumstats2 = load_sumstats(
                sumstats_file=sumstats_file2,
                chr_num=chr_num,
                start=start,
                end=end,
                allow_swapped_indel_alleles=allow_swapped_indel_alleles,
                sumstats_format=sumstats_file2_format,
            )

        # save calss attributes

        self.df_sumstats = df_sumstats
        self.df_sumstats2 = df_sumstats2
        self.sumstats_file_suffix = sumstats_file_suffix
        self.sumstats_file2_suffix = sumstats_file2_suffix
        self.sumstats_format = sumstats_file_format
        self.sumstats_file2_format = sumstats_file2_format
        self.n1 = n1
        self.n2 = n2
        self.type1 = type1
        self.type2 = type2
        self.sdY1 = sdY1
        self.sdY2 = sdY2
        self.chr_num = chr_num
        self.start = start
        self.end = end
        self.allow_swapped_indel_alleles = allow_swapped_indel_alleles

        ## used when run coloc with susie
        self.genotypes_file = genotypes_file
        self.sample_file = sample_file
        self.incl_samples = incl_samples
        self.ldstore_exe = ldstore_exe
        self.cache_dir = cache_dir
        self.cache_format = cache_format
        self.n_threads = n_threads
        self.memory = memory

    def sync_sumstats(self):
        """
        sync sumstats to the same index
        """

        intersection_snp = set(self.df_sumstats["snp"].tolist()) & set(
            self.df_sumstats2["snp"].tolist()
        )
        logging.info(
            f"before sync , sumstats1 have {self.df_sumstats.shape[0]} SNPs, qtl have {self.df_sumstats2.shape[0]} SNPs"
        )
        self.df_sumstats = self.df_sumstats[
            self.df_sumstats["snp"].isin(intersection_snp)
        ]
        self.df_sumstats2 = self.df_sumstats2[
            self.df_sumstats2["snp"].isin(intersection_snp)
        ]

        logging.info(
            f"after sync,sumstats1 left {self.df_sumstats.shape[0]} SNPs, qtl {self.df_sumstats2.shape[0]} SNPs"
        )

        assert self.df_sumstats.shape[0] > 0, "no intersection SNPs between sumstats1"
        assert self.df_sumstats2.shape[0] > 0, "no intersection SNPs between sumstats2"

    def finemap(
        self,
        num_causal_snps,
        out_dir,
        gene_col="gene_id",
        ld_file=None,
        verbose=False,
        allow_missing=False,
        save_all=False,
    ):
        """
        run coloc
        """
        if num_causal_snps > 1:
            logging.info(
                f"run coloc with {num_causal_snps} causal SNPs, this will run coloc.susie"
            )
            self.set_locus(locus_start=self.start, locus_end=self.end)
            if ld_file is None:
                ld_arr, df_ld_snps = self.get_ld_data(
                    self.start, self.end, need_bcor=False, verbose=verbose
                )
            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)

            self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)

            # not supported for hess currently

        self.df_sumstats.rename(
            columns={
                "SNP": "snp",
                "CHR": "chrom",
                "BP": "pos",
                "A1": "alt",
                "A2": "ref",
                "P": "pvalue",
            },
            inplace=True,
        )  # polyfun format => pheweb
        # self.df_sumstats.reset_index(drop=True, inplace=True)
        # turn self.df_sumstats from polyfun format => coloc format
        self.df_sumstats = sumstats2coloc(
            self.df_sumstats, map_dict=gwas_coloc_map_dict["pheweb"]
        )
        if self.sumstats_file2_format == "pheweb":  # polyfun format => pheweb
            self.df_sumstats2.rename(
                columns={
                    "SNP": "snp",
                    "CHR": "chrom",
                    "BP": "pos",
                    "A1": "alt",
                    "A2": "ref",
                    "P": "pvalue",
                },
                inplace=True,
            )  # polyfun format => pheweb
            # self.df_sumstats2.reset_index(drop=True, inplace=True)
            self.df_sumstats2 = sumstats2coloc(
                self.df_sumstats2, map_dict=gwas_coloc_map_dict["pheweb"]
            )

        # drop duplicated SNPs
        sumstats1_duplicated_num = self.df_sumstats.duplicated(subset=["snp"]).sum()

        sumstats2_drop_duplicated_cols = ["snp"]
        if gene_col is not None:
            if gene_col in self.df_sumstats2.columns:
                sumstats2_drop_duplicated_cols.append(gene_col)

        sumstats2_duplicated_num = self.df_sumstats2.duplicated(
            subset=sumstats2_drop_duplicated_cols
        ).sum()
        if sumstats1_duplicated_num > 0:
            logging.warning(
                f"sumstats1 have {sumstats1_duplicated_num} duplicated SNPs and will be dropped"
            )
        if sumstats2_duplicated_num > 0:
            logging.warning(
                f"sumstats2 have {sumstats2_duplicated_num} duplicated SNPs and will be dropped"
            )

        self.df_sumstats = self.df_sumstats.drop_duplicates(subset=["snp"])
        self.df_sumstats2 = self.df_sumstats2.drop_duplicates(
            subset=sumstats2_drop_duplicated_cols
        )

        assert (
            self.df_sumstats.shape[0] > 0 and self.df_sumstats2.shape[0] > 0
        ), "no valid SNPs to run coloc after drop duplicated SNPs"

        # sync sumstats
        self.sync_sumstats()

        logging.info(f"df_sumstats2 finally inlcude {self.df_sumstats2.shape[0]} rows")

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri

        pandas2ri.activate()

        r = robjects.r
        rscript_path = str(Path(__file__).parent / "coloc.R")
        r.source(rscript_path)
        # if (
        #     gene_col in self.df_sumstats2.columns
        #     and "rsid" in self.df_sumstats2.columns
        #     and "maf" in self.df_sumstats2.columns
        # ):
        #     logging.info(
        #         f"{gene_col} in sumstats2, will be used as xQTL do, groupby {gene_col} and run coloc"
        #     )
        #     logging.info(f"maf and rsid is founded in sumstats2, this is ok.")
        #     logging.info(
        #         f"run coloc with {self.df_sumstats2[gene_col].unique().shape[0]} genes, {self.df_sumstats2[gene_col].unique().tolist()[:5]} as xQTL"
        #     )

        # else:
        #     # logging.info(f"{gene_col} not in sumstats2, will be used as normal sumstats and run coloc")
        #     raise NotImplementedError(
        #         f"{gene_col} or rsid or maf not in sumstats2, currently only supported for xQTL with gene_id col in it "
        #     )

        Path(out_dir).mkdir(exist_ok=True, parents=True)
        assert hasattr(
            r, "coloc_analysis"
        ), f"coloc.R script is not loaded correctly, please check the file {rscript_path}"

        # run coloc
        from rpy2.robjects.conversion import localconverter

        import rpy2.robjects.conversion as cv

        def _none2null(none_obj):
            return r("NULL")

        none_converter = cv.Converter("None converter")
        none_converter.py2rpy.register(type(None), _none2null)

        with localconverter(
            robjects.default_converter + pandas2ri.converter + none_converter
        ):

            if gene_col in self.df_sumstats2.columns:
                res_gene_df = []
                res_snp_df = []  # save all snp results
                # return self.df_sumstats2
                for gene_id, sub_df in self.df_sumstats2.groupby(gene_col):
                    logging.info(f"run coloc with {gene_id} as xQTL")
                    current_output = str(Path(out_dir) / gene_id)
                    Path(current_output).mkdir(exist_ok=True, parents=True)
                    sub_df = set_snpid_index(
                        sub_df,
                        allow_swapped_indel_alleles=self.allow_swapped_indel_alleles,
                        SNP="snp",
                        CHR="chr",
                        BP="position",
                        A1="alt",
                        A2="ref",
                    )

                    interset_snps = list(
                        set(self.df_sumstats.index.tolist())
                        & set(sub_df.index.tolist())
                    )
                    if num_causal_snps > 1:
                        interset_snps = list(
                            set(interset_snps) & set(self.df_ld.index.tolist())
                        )
                    if len(interset_snps) == 0:
                        raise ValueError(
                            f"no intersection SNPs between sumstats1 and sumstats2 for gene {gene_id}"
                        )
                    current_df1 = self.df_sumstats.loc[interset_snps]
                    current_df2 = sub_df.loc[interset_snps]

                    current_ld_df = (
                        self.df_ld.loc[interset_snps, interset_snps]
                        if num_causal_snps > 1
                        else None
                    )

                    assert current_df1.shape[0] > 0 and current_df2.shape[0] > 0

                    if save_all:
                        current_df1.to_csv(
                            str(Path(current_output) / f"sumstats1.csv.gz"),
                            compression="gzip",
                        )
                        current_df2.to_csv(
                            f"{current_output}/sumstats2.csv.gz",
                            compression="gzip",
                        )
                        if num_causal_snps > 1:
                            self.df_ld.loc[interset_snps, interset_snps].to_csv(
                                str(Path(current_output) / f"ld.csv.gz"),
                                compression="gzip",
                            )

                    if num_causal_snps > 1:
                        assert (
                            current_ld_df.shape[0]
                            == current_ld_df.shape[1]
                            == current_df1.shape[0]
                            == current_df2.shape[0]
                        )
                        assert (
                            current_ld_df.index.tolist()
                            == current_ld_df.columns.tolist()
                            == current_df1.index.tolist()
                            == current_df2.index.tolist()
                        )
                    logging.info(f"run coloc with {current_df1.shape[0]} SNPs")
                    current_res = r.coloc_analysis(
                        sumstats1=current_df1,
                        sumstats2=current_df2,
                        output_folder=current_output,
                        n1=self.n1,
                        n2=self.n2,
                        sdY1=self.sdY1,
                        sdY2=self.sdY2,
                        type1=self.type1,
                        type2=self.type2,
                        ld_df=current_ld_df,
                        num_causal_snps=num_causal_snps,
                        sumstats1_suffix=self.sumstats_file_suffix,
                        sumstats2_suffix=self.sumstats_file2_suffix,
                    )

                    # print(current_res)
                    # with robjects.conversion.localconverter(
                    #     robjects.default_converter + pandas2ri.converter
                    # ):
                    print(current_res)
                    print(type(current_res))
                    current_res_gene_df = current_res["coloc_res"]
                    current_res_gene_df["gene_id"] = gene_id
                    # print(current_res["coloc_snp_res"].column)
                    current_res_snp_df = current_res["coloc_snp_res"][
                        [
                            "snp",
                            f"varbeta_{self.sumstats_file2_suffix}",
                            f"pvalue_{self.sumstats_file2_suffix}",
                            f"beta_{self.sumstats_file2_suffix}",
                            "SNP.PP.H4",
                        ]
                    ]
                    current_res_snp_df.columns = [
                        "snp",
                        f"varbeta_{self.sumstats_file2_suffix}_{gene_id}",
                        f"pvalue_{self.sumstats_file2_suffix}_{gene_id}",
                        f"beta_{self.sumstats_file2_suffix}_{gene_id}",
                        f"SNP.PP.H4_{gene_id}",
                    ]

                    res_gene_df.append(current_res_gene_df)
                    res_snp_df.append(current_res_snp_df)

                from functools import reduce
                import pandas as pd

                res_gene_df = pd.concat(res_gene_df)
                res_snp_df = reduce(
                    lambda x, y: pd.merge(x, y, on="snp", how="outer"),
                    res_snp_df,
                )
                res_snp_df = self.df_sumstats.rename(
                    columns={
                        "beta": f"beta_{self.sumstats_file_suffix}",
                        "varbeta": f"varbeta_{self.sumstats_file_suffix}",
                        "pvalue": f"pvalue_{self.sumstats_file_suffix}",
                    }
                ).merge(res_snp_df, on="snp", how="outer")

            else:

                logging.info(
                    f"run coloc with two sumstats files: {self.sumstats_file_suffix} and {self.sumstats_file2_suffix}"
                )
                current_output = out_dir

                Path(current_output).mkdir(exist_ok=True, parents=True)

                if save_all:
                    self.df_sumstats.to_csv(
                        str(Path(current_output) / f"sumstats1.csv.gz"),
                        compression="gzip",
                    )
                    self.df_sumstats2.to_csv(
                        f"{current_output}/sumstats2.csv.gz",
                        compression="gzip",
                    )
                    if num_causal_snps > 1:
                        self.df_ld.to_csv(
                            str(Path(current_output) / f"ld.csv.gz"), compression="gzip"
                        )
                if num_causal_snps > 1:
                    assert (
                        self.df_ld.shape[0]
                        == self.df_ld.shape[1]
                        == self.df_sumstats.shape[0]
                        == self.df_sumstats2.shape[0]
                    )
                    assert (
                        self.df_ld.index.tolist()
                        == self.df_ld.columns.tolist()
                        == self.df_sumstats.index.tolist()
                        == self.df_sumstats2.index.tolist()
                    )
                res = r.coloc_analysis(
                    sumstats1=self.df_sumstats,
                    sumstats2=self.df_sumstats2,
                    output_folder=current_output,
                    n1=self.n1,
                    n2=self.n2,
                    sdY1=self.sdY1,
                    sdY2=self.sdY2,
                    type1=self.type1,
                    type2=self.type2,
                    num_causal_snps=num_causal_snps,
                    ld_df=self.ld_df if num_causal_snps > 1 else None,
                    sumstats1_suffix=self.sumstats_file_suffix,
                    sumstats2_suffix=self.sumstats_file2_suffix,
                )

                res_gene_df = res["coloc_res"]
                res_snp_df = res["coloc_snp_res"]

        res_gene_df.to_csv(Path(out_dir) / "coloc_gene.csv", index=False)
        res_snp_df.to_csv(Path(out_dir) / "coloc_snp.csv", index=False)
        return res_gene_df, res_snp_df
