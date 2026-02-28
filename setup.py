from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="BioChem",
    version="1.0.0",
    author="Arif Maulana Azis",
    author_email="titandigitalsoft@gmail.com",
    description="Library for Bioinformatics and Computational Chemistry",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Arifmaulanaazis/BioChem",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.7",
    install_requires=[
        "requests>=2.26.0",
        "beautifulsoup4>=4.10.0",
        "pandas>=1.3.0",
        "rich>=12.0.0",
        "rdkit>=2022.03.1",
        "Pillow>=9.0.0",
    ],
) 