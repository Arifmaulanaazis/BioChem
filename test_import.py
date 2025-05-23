"""
Test import untuk library BioChem
"""

import BioChem
from BioChem import AdmetLabScraper, KnapsackScraper, ProtoxScraper, MolsoftScraper

print("Import berhasil!")
print(f"Versi: {BioChem.__version__}")
print("\nKelas yang tersedia:")
print("- AdmetLabScraper")
print("- KnapsackScraper")
print("- ProtoxScraper")
print("- MolsoftScraper") 