"""
BioChem.scrapers.protox - Module for accessing toxicity data from ProTox-II

This module provides a class for retrieving and processing toxicity data of chemical compounds
from the ProTox-II website (https://tox.charite.de/protox3), a web server for
predicting chemical compound toxicity.

Main class:
- ProtoxScraper: Retrieves toxicity information of chemical compounds from ProTox-II

Example usage:
```python
from BioChem import ProtoxScraper

# Initialize scraper
scraper = ProtoxScraper(max_workers=4, auto_resume=True, wait_minutes=10)

# Retrieve data for a single SMILES
result = scraper.run("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(result)

# Retrieve data for multiple SMILES at once
smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CCO", "C1CCCCC1"]
results = scraper.run(smiles_list)
print(results)
```
"""

import requests
from bs4 import BeautifulSoup
import pandas as pd
from typing import Union, List, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import rdmolfiles

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.logging import RichHandler
import logging
import time

# Configure rich logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console, show_time=True, markup=True)],
)
logger = logging.getLogger("protox_scraper")


class ProtoxScraper:
    """
    Class for retrieving toxicity data from the ProTox-II website.
    
    This class provides functionality to submit chemical compounds in SMILES format
    to the ProTox-II server and get toxicity prediction results in DataFrame format.
    
    Attributes:
        BASE_URL (str): Base URL of ProTox-II for similarity-based search.
        max_workers (int): Maximum number of thread workers for parallel processing.
        auto_resume (bool): Whether to automatically resume after a rate limit.
        wait_minutes (int): How many minutes to wait after a rate limit.
    """
    
    BASE_URL = "https://tox.charite.de/protox3/index.php?site=compound_search_similarity"

    def __init__(self, max_workers: int = 4, auto_resume: bool = False, wait_minutes: int = 10):
        """
        Initialize ProtoxScraper.
        
        Args:
            max_workers (int, optional): Maximum number of thread workers. Default 4.
            auto_resume (bool, optional): Whether to automatically resume after a rate limit. Default False.
            wait_minutes (int, optional): How many minutes to wait after a rate limit. Default 10.
        """
        self.max_workers = max_workers
        self.auto_resume = auto_resume
        self.wait_minutes = wait_minutes


    def smiles_to_molblock(self, smiles: str) -> str:
        """
        Convert SMILES to MolBlock format.
        
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
        return rdmolfiles.MolToMolBlock(mol)

    def fetch_html(self, smiles: str) -> str:
        """
        Fetch HTML data from ProTox-II for the given SMILES.
        
        Args:
            smiles (str): Valid SMILES string.
            
        Returns:
            str: HTML content of the request result.
            
        Raises:
            ConnectionError: If connection to server fails.
            RuntimeError: If rate limited and auto_resume is False.
        """
        molblock = self.smiles_to_molblock(smiles)
        payload = {
            "smilesString": molblock,
            "defaultName": "Tamoxifen",
            "smiles": smiles,
            "pubchem_name": ""
        }

        response = requests.post(self.BASE_URL, data=payload, timeout=30)
        if response.status_code != 200:
            raise ConnectionError(f"Failed to retrieve data for {smiles}, status: {response.status_code}")

        # Check if blocked due to rate limit
        if "You reached the limit of allowed queries" in response.text:
            logger.warning(f"[bold yellow]⚠ Hit rate limit while requesting: {smiles}[/]")

            if self.auto_resume:
                logger.info(f"[blue]⏳ Waiting {self.wait_minutes} minutes before continuing...[/]")
                time.sleep(self.wait_minutes * 60)
                return self.fetch_html(smiles)  # Retry after sleeping
            else:
                raise RuntimeError("Rate limit detected. Try again later.")

        return response.text

    def parse_html(self, html: str, smiles: str) -> pd.DataFrame:
        """
        Parse HTML from ProTox-II results.
        
        Args:
            html (str): HTML content of the request result.
            smiles (str): SMILES used in the request.
            
        Returns:
            pandas.DataFrame: DataFrame containing toxicity prediction results.
        """
        soup = BeautifulSoup(html, "html.parser")

        def extract_text(label: str) -> str:
            h1 = soup.find("h1", string=lambda s: s and label in s)
            return h1.get_text(strip=True).split(":")[-1].strip() if h1 else None

        data = {
            "SMILES": [smiles],
            "Predicted LD50": [extract_text("Predicted LD50")],
            "Toxicity Class": [extract_text("Predicted Toxicity Class")],
            "Average Similarity": [extract_text("Average similarity")],
            "Prediction Accuracy": [extract_text("Prediction accuracy")]
        }

        return pd.DataFrame(data)

    def process_single(self, smiles: str) -> pd.DataFrame:
        """
        Process a single SMILES and retrieve its toxicity data.
        
        Args:
            smiles (str): Valid SMILES string.
            
        Returns:
            pandas.DataFrame: DataFrame containing toxicity prediction results.
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
        Run toxicity prediction for the given SMILES.
        
        Args:
            smiles_input (str or list): Single SMILES or list of SMILES.
            
        Returns:
            pandas.DataFrame: DataFrame containing toxicity prediction results for all SMILES.
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

                    # Update progress bar
                    progress.update(task_id, advance=1, description=f"[yellow]⏳ Processing SMILES... [{completed}/{total}]")

        final_df = pd.concat(results, ignore_index=True)
        logger.info(f"[bold green]🏁 Finished![/] Total successful: {len(final_df)} / {len(smiles_list)}")
        return final_df 