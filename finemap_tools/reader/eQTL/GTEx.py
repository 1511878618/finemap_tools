"""
eQTL of GTEx v8 can be found:[eQTL catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)

"""

from finemap_tools.reader import tabix_reader
import pandas as pd

def GTEx_tabix_reader(tabix_file, region):
    """
    Read data from a tabix file for a specific region.

    Parameters:
    tabix_file (str): The path to the tabix file.
    region (str): The region to extract data from in the format "chromosome:start-end".

    Returns:
    pandas.DataFrame or None: The extracted data as a pandas DataFrame, or None if the data is not available.

    Example:
    data = GTEx_tabix_reader(tabix_file="Kidney_Cortex.tsv.gz", region="1:15573000-16174000")
    """
    data = tabix_reader(file_path=tabix_file, region=region)
    if data is None:
        return None
    colnames = pd.read_csv(tabix_file, sep="\t", nrows=0).columns
    data.columns = colnames

    return data