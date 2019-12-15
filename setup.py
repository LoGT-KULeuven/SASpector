import os
import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'SASpector',
    version = '0.0.1',
    author = 'KU Leuven - IBP Group',
    author_email = 'alecorrojo@gmail.com',
    description = 'Short-read Assembly inSpector',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/alejocrojo09/IBP19-20',
    packages = setuptools.find_packages(),
    keywords = '',
    install_requires = [
        'progressbar',
        'pandas',
        'seaborn',
        'matplotlib',
        'biopython',
        'sourmash',
        'numpy'
    ],
    data_files = [('', ['SASpector/saspector_proteindb.fasta'])],
    scripts = ['SASpector/SASpector', 'SASpector/SASpectorcheck', 'SASpector/coverage.py', 'SASpector/gene_predict.py', 'SASpector/kmer.py', 'SASpector/quastunmap.py', 'SASpector/mapper.py', 'SASpector/summary.py'],
    classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: Unix",
    "Operating System :: MacOS :: MacOS X"
    ],
    include_package_data = True,
    zip_safe = False,
    python_requires='>=3.4'

)
