import os
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'SASpector',
    version = '0.0.5',
    author = 'KU Leuven',
    author_email = 'alecorrojo@gmail.com, cedric.lood@kuleuven.be',
    maintainer = 'Computational systems biology group',
    maintainer_email = 'vera.vannoort@kuleuven.be',
    description = 'Short-read Assembly inSpector',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/LoGT-KULeuven/SASpector',
    packages = setuptools.find_packages(),
    keywords = '',
    install_requires = [
        'progressbar',
        'pandas',
        'seaborn>=0.10',
        'matplotlib',
        'biopython',
        'sourmash',
        'numpy'
    ],
    scripts = [
        'SASpector/SASpector',
        'SASpector/coverage.py',
        'SASpector/gene_predict.py',
        'SASpector/kmer.py',
        'SASpector/quastunmap.py',
        'SASpector/mapper.py',
        'SASpector/summary.py',
        'SASpector/tandem_repeats.py',
        'SASpector/select_mash.py'
    ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X"
    ],
    zip_safe = False,
    python_requires='>=3.4'
)
