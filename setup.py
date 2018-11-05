#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
from tenxt import __version__

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['Click>=6.0',
                'numpy',
                'scipy',
                'future' ]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author='Kwangbom \"KB\" Choi, Ph.D.',
    author_email='kb.choi@jax.org',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="A versatile python framework for single-cell sequencing data analysis",
    entry_points={
        'console_scripts': [
            'tenxt=tenxt.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='tenxt',
    name='tenxt',
    packages=find_packages(include=['tenxt']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/churchill-lab/tenxt',
    version=__version__,
    zip_safe=False,
)
