
import pandas as pd 

def iter_count(file_name):
    from itertools import takewhile, repeat
    buffer = 1024 * 1024
    with open(file_name) as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        return sum(buf.count('\n') for buf in buf_gen)

def check_comment_lines(file, comment="##"):
    with open (file, "r") as f:
        idx = 0
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith(comment):
                break
            idx += 1
    return idx


def add_ID(df: pd.DataFrame, col_list: list, sort: bool = False, new_col: str = "ID", id_sep: str = ":", inplace: bool = False) -> pd.DataFrame:
    """
    Adds an ID column to a DataFrame based on specified columns.

    Args:
        df (pandas.DataFrame): The DataFrame to add the ID column to.
        col_list (list): A list of column names to use for creating the ID.
        sort (bool, optional): Whether to sort the A1 and A2 values before creating the ID. Defaults to False.
        new_col (str, optional): The name of the new ID column. Defaults to "ID".
        id_sep (str, optional): The separator to use between the ID components. Defaults to ":".
        inplace (bool, optional): Whether to modify the DataFrame in-place or return a new DataFrame. Defaults to False.

    Returns:
        pandas.DataFrame or None: If inplace is True, returns None. Otherwise, returns a new DataFrame with the ID column added.
    Examples:
        eQTL_df['ID_sorted'] = add_ID(eQTL_df, ['chromosome', 'position', 'ref', 'alt'], sort=True, inplace=False)
        gwas['ID'] = add_ID(gwas, col_list=["chrom", "pos", "ref", "alt"])

    """
    
    def func(df: pd.DataFrame) -> str:
        chr = str(df[col_list[0]])
        pos = str(df[col_list[1]])
        A1 = df[col_list[2]]
        A2 = df[col_list[3]]
        if sort: 
            sorted_A1_A2 = sorted([A1, A2])
            return chr + id_sep + pos + id_sep + sorted_A1_A2[0] + id_sep + sorted_A1_A2[1]
        else:
            return chr + id_sep + pos + id_sep + A1 + id_sep + A2
    if inplace:
        df[new_col] = df.apply(func, axis=1)
    else:
        return df.apply(func, axis=1)
        
def sort_ID(df: pd.DataFrame, id_col: str, id_sep: str = ":", inplace: bool = False) -> pd.Series:
    """
    Sorts the IDs in a DataFrame column based on chromosome, position, and alleles.

    Args:
        df (pandas.DataFrame): The DataFrame containing the IDs.
        id_col (str): The name of the column containing the IDs.
        id_sep (str, optional): The separator used in the IDs. Defaults to ":".
        inplace (bool, optional): Whether to modify the DataFrame in-place. Defaults to False.

    Returns:
        pandas.Series or None: If inplace is False, returns a new Series with sorted IDs. If inplace is True, modifies the DataFrame in-place.
    Examples:
        pvar_df['ID_sorted'] = sort_ID(pvar_df, "ID", inplace=False)
    """
    def func(id: str) -> str:
        chr, pos, A1, A2 = id.split(id_sep)
        A1_A2 = sorted([A1, A2])
        return id_sep.join([chr, pos, A1_A2[0], A1_A2[1]])
    if inplace:
        df[id_col] = df[id_col].apply(func)
    else:
        return df[id_col].apply(func)
    

# def add_ID(df:pd.DataFrame, col_list, sort=False, new_col="ID", id_sep=":", inplace=False):
#     """
#     Adds an ID column to a DataFrame based on specified columns.

#     Args:
#         df (pandas.DataFrame): The DataFrame to add the ID column to.
#         col_list (list): A list of column names to use for creating the ID.
#         sort (bool, optional): Whether to sort the A1 and A2 values before creating the ID. Defaults to False.
#         new_col (str, optional): The name of the new ID column. Defaults to "ID".
#         id_sep (str, optional): The separator to use between the ID components. Defaults to ":".
#         inplace (bool, optional): Whether to modify the DataFrame in-place or return a new DataFrame. Defaults to False.

#     Returns:
#         pandas.DataFrame or None: If inplace is True, returns None. Otherwise, returns a new DataFrame with the ID column added.
#     Examples:
#         eQTL_df['ID_sorted'] = add_ID(eQTL_df, ['chromosome', 'position', 'ref', 'alt'], sort=True, inplace=False)
#     """
    
#     def func(df):
#         chr = str(df[col_list[0]])
#         pos = str(df[col_list[1]])
#         A1 = df[col_list[2]]
#         A2 = df[col_list[3]]
#         if sort: 
#             sorted_A1_A2 = sorted([A1, A2])
#             return chr + id_sep + pos + id_sep + sorted_A1_A2[0] + id_sep + sorted_A1_A2[1]
#         else:
#             return chr + id_sep + pos + id_sep + A1 + id_sep + A2
#     if inplace:
#         df[new_col] = df.apply(func, axis=1)
#     else:
#         return df.apply(func, axis=1)
        
# def sort_ID(df, id_col, id_sep=":", inplace=False):
#     """
#     Sorts the IDs in a DataFrame column based on chromosome, position, and alleles.

#     Args:
#         df (pandas.DataFrame): The DataFrame containing the IDs.
#         id_col (str): The name of the column containing the IDs.
#         id_sep (str, optional): The separator used in the IDs. Defaults to ":".
#         inplace (bool, optional): Whether to modify the DataFrame in-place. Defaults to False.

#     Returns:
#         pandas.Series or None: If inplace is False, returns a new Series with sorted IDs. If inplace is True, modifies the DataFrame in-place.
#     Examples:
#         pvar_df['ID_sorted'] = sort_ID(pvar_df, "ID", inplace=False)
#     """

#     def func(id):
#         chr, pos, A1, A2 = id.split(id_sep)
#         A1_A2 = sorted([A1, A2])
#         return id_sep.join([chr, pos, A1_A2[0], A1_A2[1]])

#     if inplace:
#         df[id_col] = df[id_col].apply(func)
#     else:
#         return df[id_col].apply(func)
