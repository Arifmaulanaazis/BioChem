# BioChem

Integrated library for Bioinformatics and Computational Chemistry that simplifies data retrieval from various sources.

## Main Features

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
pip install BioChem-0.1.0-py3-none-any.whl
```

## Requirements

- Python 3.7+
- RDKit
- Pandas
- Beautiful Soup 4
- Requests
- Rich

## Usage

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