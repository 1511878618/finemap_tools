from typing import overload, List 
import pandas as pd


@overload 
def filter_pipline(snplist:List[str], id_sep:str, pipline:List[str])->List[str]:...
@overload
def filter_pipline(sumstats:pd.DataFrame, id_col:str, pipline:List[str])->pd.DataFrame:...

def filter_pipline(snplist:List[str]=None, id_sep:str=":", pipline:List[str]=None, sumstats:pd.DataFrame=None, id_col:str=None):
    if sumstats is not None and id_col is not None and id_sep is not None:
        return filter_pipline_by_sumstats(sumstats, id_col, id_sep, pipline)
    elif snplist is not None and id_sep is not None:
        return filter_pipline_by_snplist(snplist, id_sep, pipline)


def filter_pipline_by_sumstats(sumstats:pd.DataFrame=None, id_col:str=None, id_sep:str=":",pipline:List[str]=None) -> pd.DataFrame:
    snplist = sumstats[id_col].tolist()
    passed_snplist = filter_pipline(snplist, id_sep=":", pipline=pipline)
    return sumstats[sumstats[id_col].isin(passed_snplist)]

def filter_pipline_by_snplist(snplist, id_sep=":", pipline=None):
    """
    pipline:
    1. drop ambiguous alleles  ambiguous_alleles
    2. drop biallelic snps  biallelic
    """

    if pipline == None:
        pipline = ['ambiguous_alleles', "biallelic"]
    
    if "ambiguous_alleles" in pipline:
        ambiguous_alleles = find_amibiguous_alleles(snplist, id_sep)
        print(f"drop {len(ambiguous_alleles)} ambiguous alleles")
        snplist = list(set(snplist) - set(ambiguous_alleles))
    if "biallelic" in pipline:
        bi_allelic = find_Biallelic_snp(snplist, id_sep)
        print(f"drop {len(bi_allelic)} biallelic snps")
        snplist = list(set(snplist) - set(bi_allelic))

    return snplist

@overload
def find_amibiguous_alleles(snp_id_list:List[str], id_sep:str) -> pd.DataFrame:...  
@overload
def find_amibiguous_alleles(snplist_df:pd.DataFrame, pos_col:str, id_col:str, chr_col:str) -> pd.DataFrame:...

def find_amibiguous_alleles(snp_id_list:List[str]=None, id_sep:str=None,
                            snplist_df:pd.DataFrame=None, pos_col:str=None, id_col:str=None, chr_col:str=None):
    if snp_id_list is not None and id_sep is not None:
        snplist_df = pd.DataFrame([[i] + i.split(id_sep) for i in snp_id_list], columns=["ID","chr", "pos", "A1", "A2"])
        pos_col = "pos"
        id_col = "ID"
        chr_col = "chr"
    elif snplist_df is not None and pos_col is not None and id_col is not None:
        snplist_df = snplist_df 

    return snplist_df[snplist_df.apply(lambda x: is_ambiguous_alleles(x[id_col], x["A1"], x["A2"]), axis=1)]['ID'].tolist()

@overload
def is_ambiguous_alleles(id:str, ref:str=None, alt:str=None):...
@overload
def is_ambiguous_alleles( ref:str, alt:str,id:str=None):...

def is_ambiguous_alleles(id:str=None, ref:str=None, alt:str=None):
    """
    Check if the reference and alternate alleles are ambiguous.

    Parameters:
    - id (str): The SNP ID in the format "chr:pos:a1:a2".
    - ref (str): The reference allele.
    - alt (str): The alternate allele.

    Returns:
    - bool: True if the alleles are ambiguous, False otherwise.

    Raises:
    - ValueError: If neither id nor ref and alt are provided.
    Example:
        is_ambiguous_alleles("1:100:A:T") # False
    """
    if id is not None:
        chr, pos ,a1, a2 = id.split(":")
    elif ref is not None and alt is not None:
        a1 = ref
        a2 = alt
    else:
        raise ValueError("Please provide id or ref and alt")

    if ref == "A" and alt == "T":
        return True
    elif ref == "T" and alt == "A":
        return True
    elif ref == "C" and alt == "G":
        return True
    elif ref == "G" and alt == "C":
        return True
    else:
        return False

from typing import List

@overload 
def find_Biallelic_snp(snp_id_list:List[str], id_sep:str):...
@overload
def find_Biallelic_snp(snplist_df:pd.DataFrame, pos_col:str, id_col:str, chr_col:str):...

def find_Biallelic_snp(snp_id_list:List[str]=None,id_sep:str=None,snplist_df:pd.DataFrame=None, pos_col:str=None, id_col:str=None, chr_col:str=None):
    """
    Find biallelic SNPs from a list of SNP IDs or a DataFrame.

    Parameters:
    - snp_id_list (List[str]): A list of SNP IDs in the format "chr:pos:a1:a2".
    - id_sep (str): The separator used in the SNP IDs.
    - snplist_df (pd.DataFrame): A DataFrame containing SNP information.
    - pos_col (str): The column name for the SNP positions in the DataFrame.
    - id_col (str): The column name for the SNP IDs in the DataFrame.
    - chr_col (str): The column name for the chromosome information in the DataFrame.

    Returns:
    - pd.DataFrame: A DataFrame containing the biallelic SNPs.

    Note:
    - Either snp_id_list and id_sep or snplist_df, pos_col, id_col, and chr_col must be provided.
    Example:
        find_Biallelic_snp(snp_id_list=SNP_ID_real, id_sep=":")
        or 
        find_Biallelic_snp(snplist_df=snplist_df, pos_col="pos", id_col="ID", chr_col="chr")
    """
    if snp_id_list is not None and id_sep is not None:
        snplist_df = pd.DataFrame([[i] + i.split(id_sep) for i in snp_id_list], columns=["ID","chr", "pos", "A1", "A2"])
        pos_col = "pos"
        id_col = "ID"
        chr_col = "chr"
    elif snplist_df is not None and pos_col is not None and id_col is not None:
        snplist_df = snplist_df 

    pos_is_Biallelic = snplist_df.groupby([chr_col, pos_col]).apply(lambda x:  reset_df if (reset_df := x.drop_duplicates(subset=[id_col])).shape[0] > 1 else None).reset_index(drop=True)

    if len(pos_is_Biallelic.shape) == 0:
        return None 
    else:
        return pos_is_Biallelic['ID'].tolist()
