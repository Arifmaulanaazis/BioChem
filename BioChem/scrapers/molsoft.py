"""
BioChem.scrapers.molsoft - Module for accessing molecular property data from Molsoft

This module provides a class for retrieving and processing molecular property data
from the Molsoft website (https://www.molsoft.com/mprop/), a web tool for
predicting various important molecular properties in drug development.

Main class:
- MolsoftScraper: Retrieves and analyzes molecular property data

Example usage:
```python
from BioChem import MolsoftScraper

# Initialize scraper with 4 workers
scraper = MolsoftScraper(max_workers=4)

# Retrieve data for a single SMILES
result = scraper.run("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(result)

# Retrieve data for multiple SMILES at once
smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CCO", "C1CCCCC1"]
results = scraper.run(smiles_list)
print(results)
```
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import re
import pandas as pd
from bs4 import BeautifulSoup
from typing import Union, List
from concurrent.futures import ThreadPoolExecutor, as_completed
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.console import Console
from rich.logging import RichHandler
import logging

# Configure rich logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console, show_time=True, markup=True)],
)
logger = logging.getLogger("molsoft_scraper")


class MolsoftScraper:
    """
    Class for retrieving molecular property data from the Molsoft website.
    
    This class provides functionality to submit SMILES to Molsoft
    and get molecular property prediction results in DataFrame format.
    Properties retrieved include LogP, LogS, PSA, and other parameters
    relevant for drug development.
    
    Attributes:
        BASE_URL (str): Base URL of Molsoft for molecular properties.
        max_workers (int): Maximum number of thread workers for parallel processing.
    """
    
    BASE_URL = "https://www.molsoft.com/mprop/"

    def __init__(self, max_workers: int = 4):
        """
        Initialize MolsoftScraper.
        
        Args:
            max_workers (int, optional): Maximum number of thread workers. Default 4.
        """
        self.max_workers = max_workers

    def smiles_to_molblock(self, smiles: str) -> str:
        """
        Convert SMILES to MolBlock format with 2D coordinates.
        
        Args:
            smiles (str): Valid SMILES string.
            
        Returns:
            str: MolBlock representation of the SMILES.
            
        Raises:
            ValueError: If SMILES is invalid.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        AllChem.Compute2DCoords(mol)
        mol_block = Chem.MolToMolBlock(mol)
        return mol_block.replace("RDKit", "MOLSOFT", 1)

    def fetch_html(self, smiles: str) -> str:
        """
        Submit SMILES data to Molsoft and get HTML response.
        
        Args:
            smiles (str): Valid SMILES string.
            
        Returns:
            str: HTML content of the response.
            
        Raises:
            ConnectionError: If connection to server fails.
        """
        mol_block = self.smiles_to_molblock(smiles)
        payload = {
            "p": "",
            "sm": "",
            "jme_mol": mol_block,
            "act": "Search",
            "Calc": "Calculate Properties"
        }

        headers = {
            "Content-Type": "application/x-www-form-urlencoded"
        }

        response = requests.post(self.BASE_URL, data=payload, headers=headers)
        if response.status_code != 200:
            raise ConnectionError(f"Failed to retrieve data for {smiles}, status: {response.status_code}")
        return response.text

    def parse_html(self, html: str, smiles: str) -> pd.DataFrame:
        """
        Parse HTML from Molsoft results and extract molecular property information.
        
        Args:
            html (str): HTML content of the response.
            smiles (str): SMILES used in the request.
            
        Returns:
            pandas.DataFrame: DataFrame containing molecular properties.
        """
        soup = BeautifulSoup(html, "html.parser")
        b_tags = soup.find_all("b")

        def get_value(key):
            tag = next((b for b in b_tags if key in b.text), None)
            if tag:
                return tag.next_sibling.strip() if tag.next_sibling else None
            return None

        logs_text = get_value("MolLogS :")
        logs_val = None
        if logs_text:
            match = re.search(r"([-\d.]+)\s+\(in Log", logs_text)
            logs_val = match.group(1) if match else None

        bbb_score_text = get_value("BBB Score :")
        bbb_score = None
        if bbb_score_text:
            match = re.search(r"^\s*([-\d.]+)", bbb_score_text)
            bbb_score = match.group(1) if match else None

        data = {
            "SMILES": [smiles],
            "Molecular formula": [get_value("Molecular formula:")],
            "Molecular weight": [get_value("Molecular weight:")],
            "HBA": [get_value("Number of HBA:")],
            "HBD": [get_value("Number of HBD:")],
            "MolLogP": [get_value("MolLogP :")],
            "MolLogS": [logs_val],
            "MolPSA": [get_value("MolPSA :")],
            "MolVol": [get_value("MolVol :")],
            "pKa": [get_value("pKa of most Basic/Acidic group :")],
            "BBB Score": [bbb_score],
            "Number of stereo centers": [get_value("Number of stereo centers:")]
        }

        return pd.DataFrame(data)



    def process_single(self, smiles: str) -> pd.DataFrame:
        """
        Process a single SMILES and retrieve its molecular property data.
        
        Args:
            smiles (str): Valid SMILES string.
            
        Returns:
            pandas.DataFrame: DataFrame containing molecular properties.
        """
        try:
            html = self.fetch_html(smiles)
            df = self.parse_html(html, smiles)
            console.log(f"[green]✔ Success:[/] {smiles}")
            return df
        except Exception as e:
            console.log(f"[red]✖ Failed:[/] {smiles} | {e}")
            return pd.DataFrame()


    def run(self, smiles_input: Union[str, List[str]]) -> pd.DataFrame:
        """
        Run the process of retrieving molecular property data for given SMILES.
        
        Args:
            smiles_input (str or list): Single SMILES or list of SMILES.
            
        Returns:
            pandas.DataFrame: DataFrame containing molecular properties for all SMILES.
        """
        smiles_list = [smiles_input] if isinstance(smiles_input, str) else smiles_input
        logger.info(f"[cyan]🔍 Starting scraping {len(smiles_list)} SMILES...[/]")

        results = []
        total = len(smiles_list)
        completed = 0

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
            transient=True
        ) as progress:

            task_id = progress.add_task(f"[yellow]⏳ Processing SMILES... [0/{total}]", total=total)

            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {
                    executor.submit(self.process_single, smiles): smiles
                    for smiles in smiles_list
                }

                for future in as_completed(futures):
                    result = future.result()
                    completed += 1

                    if not result.empty:
                        results.append(result)

                    # Update progress bar & description
                    progress.update(task_id, advance=1, description=f"[yellow]⏳ Processing SMILES... [{completed}/{total}]")

        final_df = pd.concat(results, ignore_index=True)
        logger.info(f"[bold green]🏁 Finished![/] Total successful: {len(final_df)} / {len(smiles_list)}")
        return final_df