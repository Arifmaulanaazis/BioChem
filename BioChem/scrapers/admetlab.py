"""
BioChem.scrapers.admetlab - Module for accessing ADMETlab data

This module provides a class for retrieving and processing ADMET property data
(Absorption, Distribution, Metabolism, Excretion, and Toxicity) of chemical compounds
from the ADMETlab website (https://admetlab3.scbdd.com).

Main class:
- AdmetLabScraper: Retrieves and analyzes ADMET property data of chemical compounds

Example usage:
```python
from BioChem import AdmetLabScraper

# Initialize scraper with 4 workers and batch size 50
scraper = AdmetLabScraper(max_workers=4, max_batch_size=50)

# Retrieve data for a single SMILES
results = scraper.run("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(results)

# Retrieve data for multiple SMILES at once
smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CCO", "C1CCCCC1"]
results = scraper.run(smiles_list)
print(results)
```
"""

import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from io import StringIO
import pandas as pd
import logging
from rich.logging import RichHandler
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

# Setup rich logging
console = Console()

# Setup FileHandler for log file
file_handler = logging.FileHandler("logs.log", encoding="utf-8")
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(logging.Formatter(
    fmt="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="[%Y-%m-%d %H:%M:%S]"
))

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console, show_time=True, markup=True), file_handler],
)
logger = logging.getLogger("admetlab_scraper")


class AdmetLabScraper:
    """
    Class for retrieving ADMET property data from the ADMETlab website.
    
    This class provides functionality to submit SMILES to ADMETlab
    and get ADMET property prediction results in DataFrame format.
    
    Attributes:
        BASE_URL (str): Base URL of ADMETlab.
        INDEX_URL (str): URL of the main ADMETlab Screening page.
        POST_URL (str): URL endpoint for submitting SMILES requests.
        max_workers (int): Maximum number of thread workers for parallel processing.
        max_batch_size (int): Maximum number of SMILES in one batch.
    """
    
    BASE_URL = "https://admetlab3.scbdd.com"
    INDEX_URL = f"{BASE_URL}/server/screening"
    POST_URL = f"{BASE_URL}/server/screeningCal"

    def __init__(self, max_workers: int = 4, max_batch_size: int = 100):
        """
        Initialize AdmetLabScraper.
        
        Args:
            max_workers (int, optional): Maximum number of thread workers. Default 4.
            max_batch_size (int, optional): Maximum number of SMILES in one batch. Default 100.
            
        Raises:
            ValueError: If max_batch_size is not in the range 1-100.
        """
        if max_batch_size > 100:
            raise ValueError("⚠️ Maximum batch size is 100.")
        if max_batch_size < 1:
            raise ValueError("⚠️ Minimum batch size is 1.")

        self.max_workers = max_workers
        self.max_batch_size = max_batch_size


    def _get_csrf_token(self, session):
        """
        Get CSRF token from the ADMETlab page.
        
        Args:
            session (requests.Session): Active HTTP session.
            
        Returns:
            str: The obtained CSRF token.
            
        Raises:
            ValueError: If CSRF token is not found.
        """
        response = session.get(self.INDEX_URL, verify=False)
        soup = BeautifulSoup(response.text, 'html.parser')
        token_input = soup.find('input', {'name': 'csrfmiddlewaretoken'})
        if not token_input:
            raise ValueError("CSRF token not found.")
        return token_input['value']

    def _submit_smiles(self, session, smiles_text, token):
        """
        Submit SMILES to ADMETlab.
        
        Args:
            session (requests.Session): Active HTTP session.
            smiles_text (str): SMILES text to submit.
            token (str): Valid CSRF token.
            
        Returns:
            requests.Response: HTTP response from the submitted request.
        """
        headers = {
            'Referer': self.INDEX_URL,
            'Origin': self.BASE_URL,
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/135.0.0.0 Safari/537.36',
            'Content-Type': 'application/x-www-form-urlencoded',
        }
        data = {
            'csrfmiddlewaretoken': token,
            'smiles-list': smiles_text,
            'method': '2'
        }
        return session.post(self.POST_URL, headers=headers, data=data, verify=False)

    def _parse_summary(self, soup):
        """
        Parse summary of results from the ADMETlab results page.
        
        Args:
            soup (BeautifulSoup): BeautifulSoup object of the results page.
            
        Returns:
            dict: Summary of results containing molecule counts.
        """
        summary = {}
        success = 0
        invalid = 0

        # Find all card titles
        cards = soup.find_all('div', class_='info-card')

        for card in cards:
            title = card.find('h5', class_='card-title')
            if title:
                title_text = title.get_text(strip=True).lower()
                number_tag = card.find('h6')
                if number_tag:
                    number = int(number_tag.get_text(strip=True))
                    if 'success' in title_text:
                        summary['success_molecules'] = number
                        success += number
                    elif 'invalid' in title_text:
                        summary['invalid_molecules'] = number
                        invalid += number
                    elif 'total' in title_text:
                        summary['total_molecules'] = number
                        
        return summary

    def _get_csv_url(self, soup):
        """
        Get the CSV download URL from the ADMETlab results page.
        
        Args:
            soup (BeautifulSoup): BeautifulSoup object of the results page.
            
        Returns:
            str: Complete URL for downloading the CSV results file.
        """
        scripts = soup.find_all('script')
        for script in scripts:
            if script.string:
                match = re.search(r'window\.open\(["\'](.*?)\.csv["\']\)', script.string)
                if match:
                    csv_url = match.group(1) + ".csv"
                    return urljoin(self.BASE_URL, csv_url)
            

    def _process_batch(self, smiles_batch):
        """
        Process a batch of SMILES and get results.
        
        Args:
            smiles_batch (list): List of SMILES to process.
            
        Returns:
            pandas.DataFrame: DataFrame containing results for the SMILES batch.
        """
        session = requests.Session()
        try:
            token = self._get_csrf_token(session)
            smiles_text = "\r\n".join(smiles_batch)
            response = self._submit_smiles(session, smiles_text, token)
            soup = BeautifulSoup(response.text, 'html.parser')

            summary = self._parse_summary(soup)
            logger.info(f"[green]✔ Batch of {len(smiles_batch)} molecules. Invalid: {summary.get('invalid_molecules')}[/]")

            csv_url = self._get_csv_url(soup)
            csv_response = session.get(csv_url)
            df = pd.read_csv(StringIO(csv_response.text))
            return df
        except Exception as e:
            logger.error(f"[red]✖ Failed to process batch: {e}[/]")
            return pd.DataFrame()

    def run(self, smiles_input):
        """
        Run the ADMET data retrieval process for the given SMILES.
        
        Args:
            smiles_input (str or list): Single SMILES or list of SMILES.
            
        Returns:
            pandas.DataFrame: DataFrame containing ADMET property prediction results.
            
        Raises:
            TypeError: If smiles_input is not a string or list.
        """
        if isinstance(smiles_input, str):
            smiles_list = [smiles_input]
        elif isinstance(smiles_input, list):
            smiles_list = smiles_input
        else:
            raise TypeError("SMILES must be a string or list of strings.")

        chunks = [smiles_list[i:i + self.max_batch_size] for i in range(0, len(smiles_list), self.max_batch_size)]
        logger.info(f"[cyan]🔍 Starting scraping of {len(smiles_list)} SMILES in {len(chunks)} batches (max {self.max_batch_size}/batch)...[/]")

        results = []
        total = len(chunks)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
            transient=True
        ) as progress:

            task_id = progress.add_task(f"[yellow]⏳ Processing batches... [0/{total}]", total=total)

            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {executor.submit(self._process_batch, batch): i for i, batch in enumerate(chunks)}

                for future in as_completed(futures):
                    result = future.result()
                    if not result.empty:
                        results.append(result)
                    progress.update(task_id, advance=1,
                        description=f"[yellow]⏳ Processing batches... [{progress.tasks[0].completed}/{total}]")

        final_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
        logger.info(f"[bold green]🏁 Done![/] Total successful molecules: {len(final_df)}")
        return final_df 