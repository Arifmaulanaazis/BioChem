"""
Basic usage examples of the BioChem library

This script demonstrates how to use various modules from the BioChem library
to retrieve data from different sources.
"""

import pandas as pd
from BioChem import AdmetLabScraper, KnapsackScraper, ProtoxScraper, MolsoftScraper

# SMILES data for testing
smiles_list = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",  # Capsaicin
    "OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@H](O)[C@@H](O)O[C@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O",  # Sucrose
]

# Example of using ADMETlab
def test_admetlab():
    print("\n===== ADMETLAB SCRAPER =====")
    scraper = AdmetLabScraper(max_workers=2, max_batch_size=10)
    results = scraper.run(smiles_list[0])  # Only using Aspirin
    print(f"Results for Aspirin:\n{results}\n")

# Example of using KNApSAcK
def test_knapsack():
    print("\n===== KNAPSACK SCRAPER =====")
    scraper = KnapsackScraper(search_type="all", keyword="Ginkgo Biloba", max_workers=2)
    results = scraper.search()
    print(f"Results for keyword 'Ginkgo Biloba':\n{results.head(3)}\n")
    
    # Saving as an Excel file
    # results.to_excel("ginkgo_biloba_results.xlsx", index=False)

# Example of using ProTox
def test_protox():
    print("\n===== PROTOX SCRAPER =====")
    scraper = ProtoxScraper(max_workers=2)
    results = scraper.run(smiles_list[0])  # Only using Aspirin
    print(f"Toxicity results for Aspirin:\n{results}\n")

# Example of using Molsoft
def test_molsoft():
    print("\n===== MOLSOFT SCRAPER =====")
    scraper = MolsoftScraper(max_workers=2)
    results = scraper.run(smiles_list[0])  # Only using Aspirin
    print(f"Molecular properties for Aspirin:\n{results}\n")

# Example of integrating all modules
def test_integration():
    print("\n===== INTEGRATION OF ALL MODULES =====")
    
    # Using Aspirin SMILES
    aspirin = smiles_list[0]
    
    # Retrieve data from all sources
    molsoft = MolsoftScraper(max_workers=1).run(aspirin)
    protox = ProtoxScraper(max_workers=1).run(aspirin)
    admetlab = AdmetLabScraper(max_workers=1).run(aspirin)
    
    # Combine results
    print("Complete data for Aspirin from various sources:")
    print("1. Molecular properties (Molsoft):")
    print(molsoft)
    print("\n2. Toxicity (ProTox):")
    print(protox)
    print("\n3. ADMET properties (ADMETlab):")
    print(admetlab.head())
    

if __name__ == "__main__":
    print("BioChem Library Usage Examples")
    
    # Choose one function to run
    # or run all sequentially
    
    test_admetlab()
    # test_knapsack()
    # test_protox()
    # test_molsoft()
    # test_integration() 