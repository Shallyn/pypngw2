#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed, 09 Oct 2024 14:44:57 +0000

@author: Shallyn
"""
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.command.install import install
from distutils.command.clean import clean
import subprocess
import sys
import os
import shutil
from pathlib import Path

class BuildCMake(build_py):
    def run(self):
        self.__root_directory = Path(os.getcwd())
        self.__module_directory = self.__root_directory / 'pypngw2'
        build_command = self.get_finalized_command('build')
        self.__cmake_build_directory = Path(build_command.build_temp) / 'cmake'
        self.__cmake_build_directory.mkdir(parents=True, exist_ok=True)

        self.remove_stale_libraries()
        subprocess.check_call([
            'cmake',
            '-S', str(self.__module_directory),
            '-B', str(self.__cmake_build_directory),
        ])
        subprocess.check_call(['cmake', '--build', str(self.__cmake_build_directory), '--config', 'Release'])

        build_py.run(self)
        self.copy_libraries()

    def get_source_dir(self):
        return os.path.abspath(os.path.dirname(__file__))

    def get_cmake_args(self):
        cmake_args = []
        # 添加需要的 CMake 参数，例如 -DCMAKE_BUILD_TYPE=Release
        return cmake_args

    def native_library_suffix(self):
        if sys.platform == 'darwin':
            return '.dylib'
        if sys.platform == 'win32':
            return '.dll'
        return '.so'

    def library_targets(self):
        return [
            ('libbasic', 'libbasic'),
            ('libcore', 'libcore'),
            ('libprec', 'libprec'),
            ('libpnfourier', 'libpnFourier'),
            ('libpngw3', 'libpngw3'),
        ]

    def remove_stale_libraries(self):
        for subdir, _ in self.library_targets():
            for lib_file in (self.__module_directory / subdir).glob('*'):
                if lib_file.suffix in ('.so', '.dylib', '.dll', '.pyd'):
                    lib_file.unlink()

    def copy_libraries(self):
        suffix = self.native_library_suffix()
        for subdir, libname in self.library_targets():
            candidates = sorted(self.__cmake_build_directory.rglob(f'{libname}{suffix}'))
            if not candidates:
                raise FileNotFoundError(f'cannot find {libname}{suffix} under {self.__cmake_build_directory}')
            lib_file = candidates[0]

            target_dir = Path(self.build_lib) / 'pypngw2' / subdir
            target_dir.mkdir(parents=True, exist_ok=True)
            for stale_file in target_dir.glob('*'):
                if stale_file.suffix in ('.so', '.dylib', '.dll', '.pyd'):
                    stale_file.unlink()
            shutil.copy2(lib_file, target_dir / lib_file.name)

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

class InstallWithLibraries(install):
    def run(self):
        install.run(self)

setup(
    name='pypngw2',
    version='1.0.0',
    packages=find_packages(),
    cmdclass={
        'build_py': BuildCMake,
        'install': InstallWithLibraries,
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
        'pypngw2': [
            'libbasic/*.so', 'libbasic/*.dylib', 'libbasic/*.dll',
            'libcore/*.so', 'libcore/*.dylib', 'libcore/*.dll',
            'libprec/*.so', 'libprec/*.dylib', 'libprec/*.dll',
            'libpnfourier/*.so', 'libpnfourier/*.dylib', 'libpnfourier/*.dll',
            'libpngw3/*.so', 'libpngw3/*.dylib', 'libpngw3/*.dll',
        ]
    },
    zip_safe=False,
)
