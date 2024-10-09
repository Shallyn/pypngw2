#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu, 12 Sep 2024 23:22:28 +0000

@author: Shallyn
"""
from setuptools.command.build_py import build_py
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.command.clean import clean
import subprocess
import sys
import os
import shutil
import glob

class BuildCMake(build_py):
    def run(self):
        # CMake construct C code
        self.__top_directory = top_directory = os.getcwd()
        self.__build_directory = build_directory = os.path.abspath('build')
        if not os.path.exists(build_directory):
            os.makedirs(build_directory)
        # Run CMake configure
        subprocess.check_call(['cmake', '..'], cwd=build_directory)
        # subprocess.check_call(['cmake', self.get_source_dir()] + self.get_cmake_args(), cwd=self.build_temp)
        # Run make
        subprocess.check_call(['make'], cwd=top_directory)
        # subprocess.check_call(['cmake', '--build', '.', '--config', 'Release'], cwd=build_directory)
        # copy dylibs
        # self.copy_libraries()
        # run build_py
        build_py.run(self)

    def get_source_dir(self):
        return os.path.abspath(os.path.dirname(__file__))

    def get_cmake_args(self):
        cmake_args = []
        # 添加需要的 CMake 参数，例如 -DCMAKE_BUILD_TYPE=Release
        return cmake_args

    def copy_libraries(self):
        # subpaths
        libs = [
            ('libbasic', 'libbasic'),
            ('libcore', 'libcore'),
            ('libprec', 'libprec'),
        ]
        for subdir, libname in libs:
            # find dylibs
            lib_build_dir = os.path.join(self.__top_directory, subdir)
            extensions = ['.so', '.dylib', '.dll']
            found = False
            for ext in extensions:
                lib_pattern = os.path.join(lib_build_dir, libname + ext)
                lib_files = glob.glob(lib_pattern)
                print(f'{lib_pattern}: {os.path.exists(lib_pattern)}')
                if lib_files:
                    lib_file = lib_files[0]
                    found = True
                    break
            if not found:
                raise FileNotFoundError(f'cannot find {libname} ({lib_pattern}) in {lib_build_dir}')
            # target path
            dest_dir = os.path.join(os.path.abspath('.'), 'pypngw2', subdir)
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            dest_file = os.path.join(dest_dir, os.path.basename(lib_file))
            print(f'dest_file: {dest_file}')
            # copy dylibs
            #shutil.copyfile(lib_file, dest_file)

class CleanCommand(clean):
    """clean command"""
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
    version='1.0',
    packages=find_packages(),
    cmdclass={
        'build_py': BuildCMake,
        'clean': CleanCommand,
    },
    include_package_data=True,
    package_data={
        'pypngw2.libbasic': ['*'],
        'pypngw2.libcore': ['*'],
        'pypngw2.libprec': ['*'],
    },
    zip_safe=False,
)
