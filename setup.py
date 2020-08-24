#!/usr/bin/env python

import os

from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'bamtk', 'VERSION'), 'r') as f:
        return f.readline().strip()


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='bamtk',
    python_requires='>=3.6',
    version=version(),
    author='Corentin Hochart',
    author_email='corentin.hochart@uca.fr',
    maintainer='Corentin Hochart',
    maintainer_email='corentin.hochart.pro@gmail.com',
    packages=['bamtk'],
    package_data={'bamtk': ['VERSION']},
    entry_points={
        'console_scripts': [
            'bamtk = bamtk.__main__:main'
        ]
    },
    url='https://github.com/chochart/bamtk',
    license='GPL3',
    description='A toolkit for dealing with bam alignement files.',
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=["biolib>=0.1.6","numpy>=1.19.1"],
    data_files=[("", ["LICENSE"])]
)