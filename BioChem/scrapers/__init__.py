"""
BioChem.scrapers - Scraper modules for BioChem

This module contains various scraper classes for retrieving data from different web sources
relevant to biochemistry and medicinal chemistry.

Available modules:
- admetlab: Retrieves ADMET property data of chemical compounds from ADMETlab
- knapsack: Retrieves metabolite data from KNApSAcK database
- protox: Retrieves toxicity data of chemical compounds from ProTox
- molsoft: Retrieves molecular property data from Molsoft
"""

from BioChem.scrapers.admetlab import AdmetLabScraper
from BioChem.scrapers.knapsack import KnapsackScraper
from BioChem.scrapers.protox import ProtoxScraper
from BioChem.scrapers.molsoft import MolsoftScraper

__all__ = [
    'AdmetLabScraper',
    'KnapsackScraper',
    'ProtoxScraper',
    'MolsoftScraper'
] 