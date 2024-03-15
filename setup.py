# coding:utf-8
# Copyright (c) 2023  Tingfeng Xu. All Rights Reserved.

from setuptools import find_packages
from setuptools import setup

import finemap_tools
import os 
from pathlib import Path
script_path = os.path.dirname(os.path.abspath(__file__)) + "/scripts"
scripts = [str(i) for i in Path(script_path).rglob("*.py") if "finemap_tools" in str(i)]

with open("requirements.txt") as file:
    REQUIRED_PACKAGES = file.read()

setup(
    name='finemap_tools',
    version=finemap_tools.__version__.replace('-', ''),
    description=('finemap_tools'),
    long_description='',
    # url='https://github.com/1511878618/cadFace',
    author='Tingfeng Xu',
    author_email='xutingfeng@big.ac.cn',
    install_requires=REQUIRED_PACKAGES,
    packages=find_packages(),
    scripts=scripts
    )
