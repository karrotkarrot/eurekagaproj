"""
PE Visualization & Results Export Module.

Generates tabular summary DataFrames, CSV exports, and 2D molecular drawings
for evaluated PE polymer structures using RDKit and Pandas.
"""

from typing import List, Dict, Any, Optional, Tuple
import os

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError:
    Chem = None
    Draw = None

class PEVisualizer:
    """Handles visualization and reporting for PE optimization runs."""

    @staticmethod
    def create_dataframe(scored_population: List[Tuple[Any, Dict[str, Any]]]) -> Any:
        """
        Converts a list of scored PE individuals into a Pandas DataFrame or dict list.
        """
        rows = []
        for pe, props in scored_population:
            if props.get("valid", False):
                rows.append({
                    "Canonical_SMILES": props.get("canonical_smiles", ""),
                    "Mn_g_mol": props.get("mn", 0.0),
                    "Total_Carbons": props.get("total_carbons", 0),
                    "Tertiary_Carbons": props.get("num_tertiary_carbons", 0),
                    "Degradation_Risk_Index": props.get("degradation_risk_index", 0.0),
                    "Estimated_Tm_C": props.get("estimated_tm_c", 0.0),
                    "Fitness_Score": props.get("fitness_score", 0.0)
                })
        if pd is not None:
            return pd.DataFrame(rows)
        return rows

    @staticmethod
    def save_dataframe_csv(df_or_rows: Any, filepath: str) -> None:
        """Exports the results to CSV file."""
        os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
        if pd is not None and isinstance(df_or_rows, pd.DataFrame):
            df_or_rows.to_csv(filepath, index=False)
        elif isinstance(df_or_rows, list) and len(df_or_rows) > 0:
            import csv
            keys = df_or_rows[0].keys()
            with open(filepath, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=keys)
                writer.writeheader()
                writer.writerows(df_or_rows)

    @staticmethod
    def draw_top_molecules(
        scored_population: List[Tuple[Any, Dict[str, Any]]],
        output_image_path: str,
        top_k: int = 5
    ) -> bool:
        """
        Generates a grid image of the top-k PE candidate structures using RDKit.
        """
        if Chem is None or Draw is None:
            print("[Warning] RDKit is not installed; molecular grid images skipped.")
            return False
            
        mols = []
        legends = []
        
        valid_candidates = [item for item in scored_population if item[1].get("valid", False)]
        for i, (pe, props) in enumerate(valid_candidates[:top_k]):
            smiles = props.get("canonical_smiles", "")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mols.append(mol)
                legend = f"Rank {i+1} | Score: {props.get('fitness_score', 0):.3f}\nMn: {props.get('mn', 0):.1f} | Risk: {props.get('degradation_risk_index', 0):.2f}"
                legends.append(legend)
                
        if not mols:
            return False
            
        try:
            os.makedirs(os.path.dirname(os.path.abspath(output_image_path)), exist_ok=True)
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=min(3, len(mols)),
                subImgSize=(350, 300),
                legends=legends
            )
            img.save(output_image_path)
            return True
        except Exception:
            return False
