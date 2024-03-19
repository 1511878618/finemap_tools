import tempfile
from typing import List, Union, overload
import shutil 
import subprocess
import pandas as pd
from pathlib import Path
import os 
def is_plink2_installed():
    if shutil.which("plink2") is None:
        return False


@overload
def plink2_cal_freq(pgen:str, snplist:List=None, outSuffix:str=None) -> pd.DataFrame:...
@overload
def plink2_cal_freq(pgen:str, snpfile:str, outSuffix:str=None) -> pd.DataFrame:...
def plink2_cal_freq(pgen:str, snpfile:str=None,snplist:List=None, outSuffix:str=None) -> pd.DataFrame:
    """
    Calculate allele frequency using plink2.

    Parameters:
    pgen (str): The path to the pgen file.
    snplist (str): The path to the snplist file.

    Returns:
    pandas.DataFrame: The allele frequency data as a pandas DataFrame.

    Example:
        1. freq = plink2_cal_freq(pgen="data/pgen/DNAJC16_ukb", snplist=SNP_ID_real, outSuffix='test') # this will create a file named test.afreq in the current directory, change with outSuffix 
        2. use tempfile to drop the output files in a temporary directory
        with tempfile.TemporaryDirectory() as tmpdirname:
            print('created temporary directory', tmpdirname)
            freq = plink2_cal_freq(pgen="data/pgen/DNAJC16_ukb", snplist=SNP_ID_real, outSuffix=tmpdirname + '/plink2_cal_freq')
            print(freq)
    

    """
    if not is_plink2_installed():
        raise Exception("plink2 is not installed.")
    if outSuffix is None:
        outSuffix = "./plink2_cal_freq"
    outSuffix = Path(outSuffix)
    outSuffix.parent.mkdir(parents=True, exist_ok=True)

    if snplist is not None:
        snplistfile = Path(str(outSuffix) + ".snplist")
        with open(snplistfile, "w") as f:
            f.write("\n".join(snplist))
    plink2_command = f"plink2 --pfile {pgen} --freq --extract {snplistfile} --out {outSuffix} --threads 4 --memory 8496"
    result = subprocess.run(plink2_command, shell=True, check=True, capture_output=True)
    freq = pd.read_csv(outSuffix.with_suffix(".afreq"), sep="\s+")

    if snplistfile.exists():
        os.remove(snplistfile)

    return freq 


@overload
def plink2_cal_LD(
    pgen: str,
    snplist: List = None,
    outSuffix: str = None,
    thread: int = 4,
    memory: int = 10240,
) -> pd.DataFrame: ...


@overload
def plink2_cal_LD(
    pgen: str, snpfile: str, outSuffix: str = None, thread: int = 4, memory: int = 10240
) -> pd.DataFrame: ...
def plink2_cal_LD(
    pgen: str,
    snpfile: str = None,
    snplist: List = None,
    outSuffix: str = None,
    thread: int = 4,
    memory: int = 10240,
) -> pd.DataFrame:
    """
    Calculate linkage disequilibrium (LD) using PLINK2.
    this should be later version of plink2 at:PLINK v2.00a6LM 64-bit Intel (2 Mar 2024)
    Args:
        pgen (str): Path to the PLINK2 pgen file.
        snpfile (str, optional): Path to the SNP file. Defaults to None.
        snplist (List, optional): List of SNPs to extract. Defaults to None.
        outSuffix (str, optional): Output file suffix. Defaults to None.

    Returns:
        pd.DataFrame: DataFrame containing the LD matrix.

    Raises:
        Exception: If PLINK2 is not installed.
    Exmaples:
        1. LD_df = plink2_cal_LD(pgen="data/pgen/DNAJC16_ukb", snplist=final_SNP_list, outSuffix='test')

        2.
        with tempfile.TemporaryDirectory() as tmpdirname:
            print('created temporary directory', tmpdirname)
            freq = plink2_cal_LD(pgen="data/pgen/DNAJC16_ukb", snplist=SNP_ID_real, outSuffix=tmpdirname + '/test')
            print(freq)
    """

    if not is_plink2_installed():
        raise Exception("plink2 is not installed.")

    if outSuffix is None:
        outSuffix = "./plink2_cal_LD"

    outSuffix = Path(outSuffix)
    outSuffix.parent.mkdir(parents=True, exist_ok=True)

    if snplist is not None:
        snplistfile = Path(str(outSuffix) + ".snplist")
        with open(snplistfile, "w") as f:
            f.write("\n".join(snplist))
    plink2_command = f"plink2 --pfile {pgen} --extract {snplistfile} --r2-unphased square --out {outSuffix} --threads {thread} --memory {memory}"
    result = subprocess.run(plink2_command, shell=True, check=True, capture_output=True)

    ld_path = outSuffix.with_suffix(".unphased.vcor2")
    LD_martix = pd.read_csv(ld_path, sep = "\t", header = None)
    LD_variants = pd.read_csv(str(ld_path) + ".vars" , sep = "\t", header = None)
    LD_martix.columns = LD_variants[0].tolist()
    LD_martix.index = LD_variants[0].tolist()

    if snplistfile.exists():
        os.remove(snplistfile)

    return LD_martix
