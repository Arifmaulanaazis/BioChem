"""
BioChem - Library for Bioinformatics and Computational Chemistry

BioChem is an integrated library that provides various tools for
analyzing, retrieving, and processing biochemical data from various web sources.

Available modules:
- admetlab: Retrieves ADMET property data of chemical compounds from ADMETlab
- knapsack: Retrieves metabolite data from KNApSAcK database
- protox: Retrieves toxicity data of chemical compounds from ProTox
- molsoft: Retrieves molecular property data from Molsoft

Basic usage:
```python
import BioChem

# Using ADMETlab scraper
scraper = BioChem.AdmetLabScraper()
results = scraper.run(["CC(=O)OC1=CC=CC=C1C(=O)O", "CCC"])
print(results)

# Using KNApSAcK scraper
knapsack = BioChem.KnapsackScraper(search_type="all", keyword="Ginkgo Biloba")
results = knapsack.search()
print(results)
```
"""

# Import all modules
from BioChem.scrapers.admetlab import AdmetLabScraper
from BioChem.scrapers.knapsack import KnapsackScraper
from BioChem.scrapers.protox import ProtoxScraper
from BioChem.scrapers.molsoft import MolsoftScraper

# Library version
__version__ = '0.1.0'

# Module names to be exposed to users
__all__ = [
    'AdmetLabScraper', 
    'KnapsackScraper', 
    'ProtoxScraper', 
    'MolsoftScraper'
] 