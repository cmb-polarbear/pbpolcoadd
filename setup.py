# -*- coding: utf-8 -*-
# Copyright 2017 Julien Peloton
# Licensed under the GPL-3.0 License, see LICENSE file for details.

from setuptools import find_packages
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('pbpolcoadd', parent_package, top_path)
    config.add_extension('pbpolcoadd',
                         sources=['pbpolcoadd/polcoadd_f.f90'],
                         libraries=['gomp'], f2py_options=[],
                         extra_f90_compile_args=[
                             '-ffixed-line-length-1000',
                             '-O3'],
                         extra_compile_args=['-fopenmp'], extra_link_args=[''],)
    return config


if __name__ == "__main__":
    ## version
    version = '0.0.1'


    setup(
        configuration=configuration,
        version=version,
        license='GPL-3.0',
        author='Dominic Beck',
        author_email='dbeck@apc.in2p3.fr',
        description='PB1 binned map-making',
        platforms='any',
        packages=find_packages(),
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Astronomy'
        ]
    )