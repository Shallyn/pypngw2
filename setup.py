#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu, 12 Sep 2024 23:22:28 +0000

@author: Shallyn
"""
import os
import pathlib
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

from setuptools import setup, find_packages
import os
# Helper to find shared libraries
def find_shared_libs():
    libs = []
    for root, dirs, files in os.walk('.'):
        for file in files:
            if file.endswith(('.so', '.pyd', '.dylib')):
                libs.append(os.path.join(root, file))
    return libs

setup(
    name='pypngw2',
    version='0.0.1',
    description='Python wrapper for PNGW2',
    packages=find_packages(),
    include_package_data=True,
    package_data={'.': find_shared_libs()},
)


