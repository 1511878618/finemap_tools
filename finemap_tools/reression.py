import statsmodels.api as sm
from typing import List
import pandas as pd


def cal_residual(df, x:List, y:str, plus_mean=True, return_model =False):
    """
    res_df = cal_residual(a, x=['age','sex'], y='ldl_a')
    res_df
    """

    print(f"passed data have {df.shape[0]} rows")
    used_df = df.copy().dropna(subset=[y] + x)
    X = sm.add_constant(used_df[x])
    Y = used_df[y]
    print(f"used data have {used_df.shape[0]} rows after dropna")


    model = sm.OLS(Y, X).fit()

    resid = model.resid
    if plus_mean:
        resid = resid + Y.mean()
    resid.name = f"{y}_residual"

    final =  df.merge(resid, left_index=True, right_index=True, how='left')

    if return_model:
        return final, model 
    else:
        return final 
