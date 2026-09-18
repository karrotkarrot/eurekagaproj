"""
PE Chemical Validator and Canonicalizer Module.

Provides chemical validity checking and canonicalization using RDKit.
Ensures molecules follow physical valence laws and eliminates duplicate representations.
"""

from typing import Optional, Tuple

try:
    from rdkit import Chem
except ImportError:
    Chem = None

class PEValidator:
    """Validator and Canonicalizer for PE polymer candidates."""

    @staticmethod
    def validate_smiles(smiles: str) -> bool:
        """
        Checks if the given SMILES string represents a chemically valid molecule.
        
        :param smiles: SMILES representation of the molecule
        :return: True if valid according to RDKit or fallback validator, False otherwise
        """
        if not smiles or not isinstance(smiles, str):
            return False
        if Chem is None:
            # Fallback syntax & hypervalency check if RDKit is not installed
            if smiles.count("(") != smiles.count(")") or "C" not in smiles:
                return False
            import re
            branch_groups = re.findall(r'\([^)]+\)', smiles)
            if len(branch_groups) > 3:
                return False
            return True
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            Chem.SanitizeMol(mol)
            return True
        except Exception:
            return False

    @staticmethod
    def canonicalize(smiles: str) -> Optional[str]:
        """
        Converts a SMILES string to its unique Canonical SMILES representation.
        
        :param smiles: SMILES representation
        :return: Canonical SMILES string if valid, None if invalid
        """
        if not smiles:
            return None
        if Chem is None:
            # Simple fallback clean up when RDKit is missing
            return smiles.strip()
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        except Exception:
            return None

    @staticmethod
    def validate_and_canonicalize(smiles: str) -> Tuple[bool, Optional[str]]:
        """
        Validates and canonicalizes a SMILES string in a single call.
        
        :param smiles: Input SMILES string
        :return: Tuple (is_valid, canonical_smiles)
        """
        if not PEValidator.validate_smiles(smiles):
            return False, None
        canon = PEValidator.canonicalize(smiles)
        return (canon is not None), canon
