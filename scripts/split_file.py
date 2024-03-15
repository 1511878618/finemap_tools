#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
@Description:       :
@Date     :2024/03/04 11:36:41
@Author      :Tingfeng Xu
@version      :1.0
"""

try:
    import pandas as pd
except:
    import subprocess
    import sys

    print("installing pandas")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "pandas"], check=True)
        print("pandas模块安装完成。")
        import pandas as pd
    except subprocess.CalledProcessError:
        print("安装pandas模块时出错。请手动安装pandas。")
        sys.exit(1)

import argparse
import os
import textwrap
import warnings


def save_subset(df, output_folder, col, save_kwargs):
    # 分别保存不同值对应的数据到新的文件
    for value, group in df.groupby(col):
        output_file = os.path.join(output_folder, f"{value}.csv")
        group.to_csv(output_file,  **save_kwargs)


def print_value_counts(df, col):
    # 打印指定列的值计数结果
    value_counts = df[col].value_counts()
    print(value_counts)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """
            %prog is ...
            @Author: xutingfeng@big.ac.cn
            Version: 1.0
            Example code:
                ./split_file.py -i input.csv -o output_folder -c 0

                ./split_file.py -i input.csv -o output_folder -c 0 --no-header

                ./split_file.py -i input.csv -o output_folder -c 0 -l
                ...
            """
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument(
        "-o", "--output-folder", required=False, help="Output folder path", default="split_file"
    )
    parser.add_argument(
        "-c", "--column", required=True, type=int, help="Column index to use, start from 1 (not 0)"
    )
    parser.add_argument(
        "--no-header", action="store_true", help="Indicate if input file has no header"
    )
    parser.add_argument(
        "-l",
        "--value-counts",
        action="store_true",
        help="Print value counts of specified column",
    )
    parser.add_argument(
        "-t", "--sep", required=False, help="Separator of input file", default="\s+"
    ) # 默认分隔符为任意空白字符
    parser.add_argument("--out-sep", required=False, help="Separator of output file, default tab ", default="\t")
    args = parser.parse_args()

    # 创建输出文件夹
    os.makedirs(args.output_folder, exist_ok=True)

    # 读取输入文件
    if args.no_header:
        df = pd.read_csv(args.input, header=None, sep=args.sep)
    else:
        df = pd.read_csv(args.input, sep = args.sep)

    col = df.columns[args.column - 1]

    if args.value_counts:
        print_value_counts(df, col)
    else:
        save_kwargs = {"sep": args.out_sep, "index": False}
        save_subset(df, args.output_folder,col, save_kwargs)


if __name__ == "__main__":
    main()
