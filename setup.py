# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open('surpyvor/version.py').read())

setup(
    name='surpyvor',
    version=__version__,
    description='Manipulate vcf files of structural variants using SURVIVOR',
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    url='https://github.com/wdecoster/surpyvor',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore',
    packages=find_packages(),
    python_requires='>=3',
    install_requires=[],
    package_data={'surpyvor': []},
    package_dir={'surpyvor': 'surpyvor'},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'surpyvor=surpyvor.surpyvor:main',
        ],
    },
)
