# BioChem 1.0.0 Stable

Integrated library for Bioinformatics and Computational Chemistry that simplifies data retrieval from various sources and provides comprehensive cheminformatics capabilities.

## Main Features

- **Cheminformatics (New in v1.0.0)**: Physicochemical properties prediction, Lipinski's Rule of Five profiling, 3D conformer generation, MMFF94/UFF energy minimization, and 2D image generation powered by RDKit.
- **ADMETlab**: Access ADMET property data of chemical compounds from ADMETlab
- **KNApSAcK**: Access plant metabolite data from KNApSAcK database
- **ProTox-II**: Access toxicity data of chemical compounds from ProTox-II
- **Molsoft**: Access molecular property data from Molsoft

## Installation

You can install BioChem using pip:

```bash
pip install BioChem
```

Or using the wheel file:

```bash
pip install https://github.com/Arifmaulanaazis/BioChem/releases/download/v1.0.0/biochem-1.0.0-py3-none-any.whl
```

## Requirements

- Python 3.7+
- RDKit >= 2022.03.1
- Pillow >= 9.0.0
- Pandas >= 1.3.0
- Beautiful Soup 4 >= 4.10.0
- Requests >= 2.26.0
- Rich >= 12.0.0

## Usage

### Example of using ChemAnalyzer (Cheminformatics)

```python
from BioChem import ChemAnalyzer

# Initialize analyzer with a SMILES string (e.g., Aspirin)
analyzer = ChemAnalyzer(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")

# Predict physicochemical properties
props = analyzer.physicochemical_properties()
print("Properties:", props)

# Check Lipinski's Rule of Five
lipinski = analyzer.lipinski_rule_of_five()
print("Lipinski Profiling:", lipinski)

# Generate 3D conformer and minimize using MMFF94
analyzer.generate_conformer()
minimize_res = analyzer.minimize_mmff94()
print("Minimization Results:", minimize_res)

# Save the minimized 3D structure to SDF
analyzer.save_conformer("aspirin_minimized.sdf")

# Generate 2D image
analyzer.generate_2d_image("aspirin.png")
```

### Example of using ADMETlab

```python
from BioChem import AdmetLabScraper

# Initialize the scraper
scraper = AdmetLabScraper(max_workers=4, max_batch_size=50)

# Retrieve data for a single SMILES
results = scraper.run("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(results)
```

### Example of using KNApSAcK

```python
from BioChem import KnapsackScraper

# Search data based on plant name
scraper = KnapsackScraper(search_type="all", keyword="Ginkgo Biloba")
results = scraper.search()
print(results)

# Export results to Excel file
results.to_excel("Ginkgo_Biloba_Data.xlsx", index=False)
```

### Example of using ProTox-II

```python
from BioChem import ProtoxScraper

# Initialize the scraper
scraper = ProtoxScraper(max_workers=4, auto_resume=True, wait_minutes=10)

# Retrieve data for a single SMILES
result = scraper.run("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(result)
```

### Example of using Molsoft

```python
from BioChem import MolsoftScraper

# Initialize the scraper
scraper = MolsoftScraper(max_workers=4)

# Retrieve data for a single SMILES
result = scraper.run("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(result)
```

## Author

- **Arif Maulana Azis** - [Arifmaulanaazis](https://github.com/Arifmaulanaazis)

## Contact

- Email: titandigitalsoft@gmail.com
- GitHub: https://github.com/Arifmaulanaazis/BioChem

## License

MIT License
