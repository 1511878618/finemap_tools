from finemap_tools.utils import check_comment_lines
import pandas as pd

def read_pvar(pvar):
    comment_lines = check_comment_lines(pvar)
    return pd.read_csv(pvar, sep="\t", skiprows=comment_lines)

