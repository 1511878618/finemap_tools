# from .base import Fine_Mapping

from finemap_tools.reader.gwas.sumstats import load_sumstats
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
    "format": {
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
    }
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


class coloc_Wrapper(object):
    def __init__(
        self,
        sumstats_file1,
        sumstats_file1_format,
        sumstats_file2,
        sumstats_file2_format,
        n1,
        n2,
        type1,
        # type2,
        sdY1,
        # sdY2,
        chr_num,
        start,
        end,
        allow_swapped_indel_alleles=False,
    ):

        df_sumstats1 = load_summstats_coloc(
            sumstats_file1,
            map_dict=gwas_coloc_map_dict.get(sumstats_file1_format, None),
            chr_num=chr_num,
            start=start,
            end=end,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
            sumstats_format=sumstats_file1_format,
        )
        df_sumstats2 = load_summstats_coloc(
            sumstats_file2,
            map_dict=eQTL_coloc_map_dict.get(sumstats_file2_format, None),
            chr_num=chr_num,
            start=start,
            end=end,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
            sumstats_format=sumstats_file2_format,
        )

        # save calss attributes
        self.sumstats_file1 = df_sumstats1
        self.sumstats_file2 = df_sumstats2
        self.n1 = n1
        self.n2 = n2
        self.type1 = type1
        # self.type2 = type2
        self.sdY1 = sdY1
        # self.sdY2 = sdY2
        self.chr_num = chr_num
        self.start = start
        self.end = end
        self.allow_swapped_indel_alleles = allow_swapped_indel_alleles

    def sync_sumstats(self):
        """
        sync sumstats to the same index
        """

        intersection_snp = set(self.sumstats_file1["snp"].tolist()) & set(
            self.sumstats_file2["snp"].tolist()
        )
        logging.info(
            f"before sync , sumstats1 have {self.sumstats_file1.shape[0]} SNPs, qtl have {self.sumstats_file2.shape[0]} SNPs"
        )
        self.sumstats_file1 = self.sumstats_file1[
            self.sumstats_file1["snp"].isin(intersection_snp)
        ]
        self.sumstats_file2 = self.sumstats_file2[
            self.sumstats_file2["snp"].isin(intersection_snp)
        ]

        logging.info(
            f"after sync,sumstats1 left {self.sumstats_file1.shape[0]} SNPs, qtl {self.sumstats_file2.shape[0]} SNPs"
        )

        assert (
            self.sumstats_file1.shape[0] > 0
        ), "no intersection SNPs between sumstats1"
        assert (
            self.sumstats_file2.shape[0] > 0
        ), "no intersection SNPs between sumstats2"

    def finemap(self, out_dir, gene_col="gene_id"):
        """
        run coloc
        """

        # drop duplicated SNPs
        sumstats1_duplicated_num = self.sumstats_file1.duplicated(subset=["snp"]).sum()
        sumstats2_duplicated_num = self.sumstats_file2.duplicated(
            subset=["snp", gene_col]
        ).sum()
        if sumstats1_duplicated_num > 0:
            logging.warning(
                f"sumstats1 have {sumstats1_duplicated_num} duplicated SNPs and will be dropped"
            )
        if sumstats2_duplicated_num > 0:
            logging.warning(
                f"sumstats2 have {sumstats2_duplicated_num} duplicated SNPs and will be dropped"
            )

        self.sumstats_file1 = self.sumstats_file1.drop_duplicates(subset=["snp"])
        self.sumstats_file2 = self.sumstats_file2.drop_duplicates(
            subset=["snp", gene_col]
        )

        assert (
            self.sumstats_file1.shape[0] > 0 and self.sumstats_file2.shape[0] > 0
        ), "no valid SNPs to run coloc after drop duplicated SNPs"

        # sync sumstats
        self.sync_sumstats()

        logging.info(
            f"sumstats_file2 finally inlcude {self.sumstats_file2.shape[0]} rows"
        )

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri

        pandas2ri.activate()

        r = robjects.r
        rscript_path = str(Path(__file__).parent / "coloc.R")
        r.source(rscript_path)
        if (
            gene_col in self.sumstats_file2.columns
            and "rsid" in self.sumstats_file2.columns
            and "maf" in self.sumstats_file2.columns
        ):
            logging.info(
                f"{gene_col} in sumstats2, will be used as xQTL do, groupby {gene_col} and run coloc"
            )
            logging.info(f"maf and rsid is founded in sumstats2, this is ok.")
            logging.info(
                f"run coloc with {self.sumstats_file2[gene_col].unique().shape[0]} genes, {self.sumstats_file2[gene_col].unique().tolist()[:5]} as xQTL"
            )

        else:
            # logging.info(f"{gene_col} not in sumstats2, will be used as normal sumstats and run coloc")
            raise NotImplementedError(
                f"{gene_col} or rsid or maf not in sumstats2, currently only supported for xQTL with gene_id col in it "
            )
        Path(out_dir).mkdir(exist_ok=True, parents=True)
        assert hasattr(
            r, "runColocAnalysis"
        ), f"coloc.R script is not loaded correctly, please check the file {rscript_path}"

        res = r.runColocAnalysis(
            data_gwas=self.sumstats_file1,
            data_eQTLs=self.sumstats_file2,
            output_folder=out_dir,
            n_gwas=self.n1,
            n_eQTL=self.n2,
            gwas_sdY=self.sdY1,
            gwas_type="quant",
        )
        from rpy2.robjects.conversion import localconverter

        with localconverter(robjects.default_converter + pandas2ri.converter):
            res_gene_df = res[0]
            res_1 = []
            for k in res[1].keys():
                pp_h4 = res[1][k].columns[-1]
                symbol_name = pp_h4.split("_")[0]
                current_df = res[1][k][
                    ["snp", "varbeta_eQTL", "pvalue_eQTL", pp_h4]
                ].rename(
                    columns={
                        "varbeta_eQTL": f"varbeta_eQTL_{symbol_name}",
                        "pvalue_eQTL": f"pvalue_eQTL_{symbol_name}",
                    }
                )
                res_1.append(current_df)
        from functools import reduce
        import pandas as pd

        res_snp_df = reduce(
            lambda x, y: pd.merge(x, y, on="snp", how="outer"),
            res_1,
        )
        res_snp_df = self.sumstats_file1.rename(
            columns={"beta": "beta_GWAS", "varbeta": "varbeta_GWAS"}
        ).merge(res_snp_df, left_index=True, right_index=True, how="inner")

        res_gene_df.to_csv(Path(out_dir) / "coloc_gene.csv", index=False)
        res_snp_df.to_csv(Path(out_dir) / "coloc_snp.csv", index=False)
        return res_gene_df, res_snp_df
