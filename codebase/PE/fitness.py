"""
PE Fitness Testing & Property Assessment Module.

Evaluates PE candidate structures based on the dependent & independent variables:
1. Number Average Molecular Weight (Mn)
2. Branch density & tertiary carbon scission susceptibility (Thermal degradation proxy)
3. Estimated Melting Point / Thermal Stability Index
"""

from typing import Dict, Any, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    Chem = None
    Descriptors = None
    rdMolDescriptors = None

from codebase.PE.validator import PEValidator

class PEFitnessEvaluator:
    """Evaluates physical properties and fitness scores for PE candidates."""

    def __init__(self, target_mn: float = 1000.0, mn_weight: float = 0.5, stability_weight: float = 0.5):
        """
        :param target_mn: Target Number Average Molecular Weight (g/mol) to aim for or exceed.
        :param mn_weight: Weight of Molecular Weight in fitness calculation.
        :param stability_weight: Weight of Thermal Stability / Low Degradation Index.
        """
        self.target_mn = target_mn
        self.mn_weight = mn_weight
        self.stability_weight = stability_weight

    def evaluate_properties(self, smiles: str) -> Dict[str, Any]:
        """
        Calculates molecular properties for a given PE SMILES string.
        
        :param smiles: Candidate PE SMILES string
        :return: Dictionary containing calculated physical properties
        """
        is_valid, canon_smiles = PEValidator.validate_and_canonicalize(smiles)
        if not is_valid or canon_smiles is None:
            return {
                "valid": False,
                "mn": 0.0,
                "num_tertiary_carbons": 0,
                "degradation_risk_index": 1.0,
                "estimated_tm_c": 0.0,
                "fitness_score": -100.0
            }

        # Fallback calculation if RDKit is not installed
        if Chem is None or Descriptors is None:
            total_c = smiles.count("C") + smiles.count("c")
            # Molecular weight approximation for alkane: C_n H_{2n+2}
            mw = (total_c * 12.011) + ((total_c * 2 + 2) * 1.008)
            num_tertiary = smiles.count("(")
            tertiary_ratio = (num_tertiary / total_c) if total_c > 0 else 0.0
            degradation_risk_index = min(1.0, tertiary_ratio * 3.0)
            estimated_tm = max(80.0, min(145.0, 135.0 - (tertiary_ratio * 60.0) + min(15.0, (mw / 100.0))))
            mn_score = min(1.0, mw / self.target_mn)
            stability_score = (1.0 - degradation_risk_index)
            composite_score = (self.mn_weight * mn_score) + (self.stability_weight * stability_score)
            return {
                "valid": True,
                "canonical_smiles": canon_smiles,
                "mn": float(mw),
                "total_carbons": total_c,
                "num_tertiary_carbons": num_tertiary,
                "degradation_risk_index": float(degradation_risk_index),
                "estimated_tm_c": float(estimated_tm),
                "fitness_score": float(composite_score)
            }

        mol = Chem.MolFromSmiles(canon_smiles)
        
        # 1. Molecular Weight (proxy for Number Average Molecular Weight Mn)
        mw = Descriptors.MolWt(mol)
        
        # 2. Count Tertiary Carbons (Branch connection points - susceptible to thermal degradation)
        tertiary_pattern = Chem.MolFromSmarts("[CH1]([#6])([#6])[#6]")
        num_tertiary = len(mol.GetSubstructMatches(tertiary_pattern)) if tertiary_pattern else 0

        # 3. Calculate Thermal Degradation Risk Index
        total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        tertiary_ratio = (num_tertiary / total_carbons) if total_carbons > 0 else 0.0
        degradation_risk_index = min(1.0, tertiary_ratio * 3.0)

        # 4. Estimated Melting Temperature Tm (°C)
        base_tm = 135.0 - (tertiary_ratio * 60.0) + min(15.0, (mw / 100.0))
        estimated_tm = max(80.0, min(145.0, base_tm))

        # 5. Composite Fitness Score Calculation
        mn_score = min(1.0, mw / self.target_mn)
        stability_score = (1.0 - degradation_risk_index)
        
        composite_score = (self.mn_weight * mn_score) + (self.stability_weight * stability_score)

        return {
            "valid": True,
            "canonical_smiles": canon_smiles,
            "mn": float(mw),
            "total_carbons": total_carbons,
            "num_tertiary_carbons": num_tertiary,
            "degradation_risk_index": float(degradation_risk_index),
            "estimated_tm_c": float(estimated_tm),
            "fitness_score": float(composite_score)
        }

    def score(self, smiles: str) -> float:
        """Helper to return only the scalar fitness score."""
        res = self.evaluate_properties(smiles)
        return res["fitness_score"]
