from finemap_tools.reader import tabix_reader
import pandas as pd
def load_GWASFormated_file(file_path:str, region:str) -> pd.DataFrame:
    """
    Load a GWAS-formatted file.

    Parameters:
    file_path (str): The path to the GWAS-formatted file.
    region (str): The genomic region to load from the file.

    Returns:
    pandas.DataFrame: The loaded data as a pandas DataFrame.

    Example:
    data = load_GWASFormated_file("T1Mapping_Cortex_20240129.csv_firstorder_Median_all_2023_GRCh38_unionKidneys.tsv.gz", region="1:15573000-16174000")
    """
    data = tabix_reader(file_path=file_path, region=region)
    if data is None:
        return None
    # As GWASFormated file may have a tabix index with header, so 
    if isinstance(data.iloc[0, 0], str):
        data.columns = data.iloc[0]
        data = data[1:]
        print("tabix have a header, so will take the first line as header and remove it.")
  
    return data

