#!/usr/bin/env python

from setuptools import setup, find_packages

setup(

    name = 'xenophile',
    version = '1.0.0',
    author = 'Michael Leeming',
    author_email = 'm.leeming@student.unimelb.edu.au',
    url = 'https://github.com/mgleeming/Xenophile',
    license = 'LICENSE.txt',
    description = 'xenophile: Non-targeted identification of xenobiotic-protein adducts',

    packages = find_packages(),
    include_package_data=True,

    entry_points = {
        'gui_scripts' : ['xenophile = xenophile.xenophile:main']
    },

    install_requires = [
        'numpy',
        'pymzml',
        'rtree',
        'pyteomics',
        'pyqtgraph'
    ],
    dependency_links = ['http://pyqt.sourceforge.net/Docs/PyQt4/installation.html'],
    zip_safe = False
)
