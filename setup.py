#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed, 09 Oct 2024 14:44:57 +0000

@author: Shallyn
"""
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from distutils.command.clean import clean
import subprocess
import sys
import os
import shutil
import glob
from pathlib import Path

class BuildCMake(build_py):
    def run(self):
        # run CMake construct C lib
        self.__root_directory = Path(os.getcwd())
        self.__module_directory = self.__root_directory / 'pypngw2'
        # self.__build_directory = self.__module_directory / 'build'
        # if not self.__build_directory.exists():
        #     self.__build_directory.mkdir(parents=True)
        # run CMake configure
        subprocess.check_call(['cmake', '.'], cwd=self.__module_directory)
        # make
        subprocess.check_call(['make'], cwd=self.__module_directory)
        # self.copy_libraries(build_directory)
        build_py.run(self)

class CleanCommand(clean):
    """clean"""
    def run(self):
        if os.path.exists('build'):
            shutil.rmtree('build')
        if os.path.exists('dist'):
            shutil.rmtree('dist')
        for root, dirs, files in os.walk('.'):
            for file in files:
                if file.endswith(('.so', '.dylib', '.dll', '.pyd')):
                    os.remove(os.path.join(root, file))
        clean.run(self)

setup(
    name='pypngw2',
    version='1.0.0',
    packages=find_packages(),
    cmdclass={
        'build_py': BuildCMake,
        'clean':CleanCommand
        },
    install_requires=[],
    author='Xiaolin Liu',
    author_email='shallyn.liu@foxmail.com',
    description='GW generator',
    url='https://github.com/Shallyn/pypngw2.git',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License', 
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    include_package_data=True,
    package_data={
        'pypngw2.libbasic': ['*'],
        'pypngw2.libcore': ['*'],
        'pypngw2.libprec': ['*'],
    },
    zip_safe=False,
)
