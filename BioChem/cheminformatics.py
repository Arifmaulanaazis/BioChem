import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem, Draw
from typing import List, Dict, Union, Any

class ChemAnalyzer:
    """
    A class for comprehensive cheminformatics operations using RDKit.
    Provides methods for calculating physicochemical properties, Lipinski's profiling,
    conformer generation, energy minimization, and 2D visualization.
    """
    
    def __init__(self, smiles: str = None, mol: Chem.Mol = None):
        """
        Initialize the ChemAnalyzer with either a SMILES string or an RDKit Mol object.
        """
        if smiles:
            self.mol = Chem.MolFromSmiles(smiles)
            if self.mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
        elif mol is not None:
            self.mol = mol
        else:
            self.mol = None

    def load_smiles(self, smiles: str):
        """Load a molecule from a SMILES string."""
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

    def physicochemical_properties(self) -> Dict[str, Union[float, int]]:
        """Calculate basic physicochemical properties."""
        if not self.mol:
            raise ValueError("No molecule loaded.")
        
        return {
            "MolecularWeight": Descriptors.MolWt(self.mol),
            "LogP": Descriptors.MolLogP(self.mol),
            "TPSA": Descriptors.TPSA(self.mol),
            "NumHDonors": Lipinski.NumHDonors(self.mol),
            "NumHAcceptors": Lipinski.NumHAcceptors(self.mol),
            "NumRotatableBonds": Lipinski.NumRotatableBonds(self.mol)
        }

    def lipinski_rule_of_five(self) -> Dict[str, Any]:
        """
        Evaluate Lipinski's Rule of Five.
        Criteria: MW <= 500, LogP <= 5, H-Donors <= 5, H-Acceptors <= 10.
        """
        if not self.mol:
            raise ValueError("No molecule loaded.")
            
        props = self.physicochemical_properties()
        
        violations = 0
        violation_details = []
        
        if props["MolecularWeight"] > 500:
            violations += 1
            violation_details.append("MW > 500")
        if props["LogP"] > 5:
            violations += 1
            violation_details.append("LogP > 5")
        if props["NumHDonors"] > 5:
            violations += 1
            violation_details.append("H-Donors > 5")
        if props["NumHAcceptors"] > 10:
            violations += 1
            violation_details.append("H-Acceptors > 10")
            
        return {
            "violations": violations,
            "details": violation_details,
            "conclusion": "Pass" if violations <= 1 else "Fail"
        }

    def generate_conformer(self, num_confs: int = 1, random_seed: int = 42) -> List[int]:
        """Generate 3D conformers for the loaded molecule."""
        if not self.mol:
            raise ValueError("No molecule loaded.")
        
        # Add hydrogens first for accurate 3D geometry
        self.mol = Chem.AddHs(self.mol)
        conformer_ids = AllChem.EmbedMultipleConfs(
            self.mol, 
            numConfs=num_confs, 
            randomSeed=random_seed
        )
        return list(conformer_ids)

    def minimize_mmff94(self, max_iters: int = 200) -> List[Dict[str, Any]]:
        """Perform MMFF94 minimization on all generated conformers."""
        if not self.mol or self.mol.GetNumConformers() == 0:
            self.generate_conformer()
            
        results = []
        res = AllChem.MMFFOptimizeMoleculeConfs(self.mol, maxIters=max_iters, mmffVariant='MMFF94')
        for conf_id, (not_converged, energy) in enumerate(res):
            results.append({
                "conformer_id": conf_id,
                "converged": not_converged == 0,
                "energy": energy
            })
        return results

    def minimize_uff(self, max_iters: int = 200) -> List[Dict[str, Any]]:
        """Perform UFF minimization on all generated conformers."""
        if not self.mol or self.mol.GetNumConformers() == 0:
            self.generate_conformer()
            
        results = []
        res = AllChem.UFFOptimizeMoleculeConfs(self.mol, maxIters=max_iters)
        for conf_id, (not_converged, energy) in enumerate(res):
            results.append({
                "conformer_id": conf_id,
                "converged": not_converged == 0,
                "energy": energy
            })
        return results

    def save_conformer(self, filepath: str, file_format: str = "sdf", conf_id: int = -1):
        """Save the minimized (or generated) molecule to a file."""
        if not self.mol:
            raise ValueError("No molecule loaded.")
            
        if file_format.lower() == "sdf":
            writer = Chem.SDWriter(filepath)
            if conf_id == -1:
                if self.mol.GetNumConformers() > 1:
                    for cid in range(self.mol.GetNumConformers()):
                        writer.write(self.mol, confId=cid)
                else:
                    writer.write(self.mol)
            else:
                writer.write(self.mol, confId=conf_id)
            writer.close()
        elif file_format.lower() == "pdb":
            Chem.MolToPDBFile(self.mol, filepath, confId=conf_id)
        else:
            raise ValueError(f"Unsupported format: {file_format}")

    def generate_2d_image(self, filepath: str, size: tuple = (300, 300)):
        """Generate a 2D drawing of the molecule and save it as an image."""
        if not self.mol:
            raise ValueError("No molecule loaded.")
        
        # Remove hydrogens for clearer 2D drawing if they were explicitly added
        display_mol = Chem.RemoveHs(self.mol)
        
        # Ensure 2D coordinates are present
        AllChem.Compute2DCoords(display_mol)
        
        Draw.MolToFile(display_mol, filepath, size=size)

    # --- Batch Operations ---
    
    @staticmethod
    def batch_predict_properties(smiles_list: List[str]) -> List[Dict[str, Any]]:
        """Batch compute physicochemical properties and Lipinski profiling."""
        results = []
        for smi in smiles_list:
            try:
                analyzer = ChemAnalyzer(smiles=smi)
                props = analyzer.physicochemical_properties()
                lipinski = analyzer.lipinski_rule_of_five()
                results.append({
                    "smiles": smi,
                    "properties": props,
                    "lipinski": lipinski,
                    "error": None
                })
            except Exception as e:
                results.append({"smiles": smi, "error": str(e)})
        return results

    @staticmethod
    def batch_minimize(smiles_list: List[str], method: str = "mmff94", save_dir: str = None) -> List[Dict[str, Any]]:
        """Batch perform conformer generation, minimization, and optionally save SDF files."""
        results = []
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
        for i, smi in enumerate(smiles_list):
            try:
                analyzer = ChemAnalyzer(smiles=smi)
                analyzer.generate_conformer(num_confs=1)
                
                if method.lower() == "mmff94":
                    min_res = analyzer.minimize_mmff94()
                elif method.lower() == "uff":
                    min_res = analyzer.minimize_uff()
                else:
                    raise ValueError(f"Unknown minimization method: {method}")
                    
                output = {
                    "smiles": smi,
                    "minimization_results": min_res,
                    "error": None
                }
                
                if save_dir:
                    filepath = os.path.join(save_dir, f"mol_{i+1}.sdf")
                    analyzer.save_conformer(filepath, file_format="sdf")
                    output["saved_file"] = filepath
                    
                results.append(output)
            except Exception as e:
                results.append({"smiles": smi, "error": str(e)})
                
        return results
