"""
PE (Polyethylene) Polymer Representation Engine.

Models PE structural parameters according to the experimental variables:
1. Polymer chain average length (Degree of Polymerization, DP)
2. Number of chain branches (N_branch) and branch length (L_branch)
3. Branch bending/directionality (Stereochemistry / branching geometry)
"""

from typing import List, Dict, Any, Optional

try:
    import selfies as sf
except ImportError:
    sf = None

class PEBranch:
    """Represents a single alkyl branch on the polyethylene main chain."""
    def __init__(self, position: int, length: int, bend_direction: str = "right"):
        """
        :param position: Carbon index along the main chain (1-indexed)
        :param length: Length of the alkyl branch in carbon atoms (1=Methyl, 2=Ethyl, 4=Butyl, etc.)
        :param bend_direction: Directionality/stereocenter configuration ('right', 'left', 'cis', 'trans')
        """
        self.position = max(1, position)
        self.length = max(1, length)
        self.bend_direction = bend_direction

    def to_smiles_substituent(self) -> str:
        """Generates the SMILES string for this side branch."""
        # Simple alkyl branch: (C) for methyl, (CC) for ethyl, (CCCC) for butyl, etc.
        # Stereochemistry markers (@/@@) can be added based on bend_direction
        ch_str = "C" * self.length
        if self.bend_direction == "left":
            return f"([C@H]{ch_str[1:] if self.length > 1 else ''})"
        elif self.bend_direction == "right":
            return f"([C@@H]{ch_str[1:] if self.length > 1 else ''})"
        else:
            return f"({ch_str})"


class PERepresentation:
    """
    Polyethylene Structure Container.
    Maintains backbone length and a list of PEBranch objects.
    """
    def __init__(self, dp: int = 10, branches: Optional[List[PEBranch]] = None):
        """
        :param dp: Degree of Polymerization (number of -CH2-CH2- repeat units in backbone).
                   Total backbone carbons = dp * 2.
        :param branches: List of PEBranch objects attached to the main chain.
        """
        self.dp = max(2, dp)
        self.total_carbons = self.dp * 2
        self.branches = branches if branches is not None else []
        self._sort_branches()

    def _sort_branches(self):
        """Ensure branches are within bounds and sorted by backbone position."""
        valid_branches = []
        for b in self.branches:
            # Keep branch position within 2 and total_carbons - 1
            pos = max(2, min(self.total_carbons - 1, b.position))
            valid_branches.append(PEBranch(pos, b.length, b.bend_direction))
        self.branches = sorted(valid_branches, key=lambda x: x.position)

    def to_smiles(self) -> str:
        """
        Constructs a linear/branched PE SMILES string by assembling backbone carbons
        and inserting branch substituents at designated carbon positions.
        """
        # Create carbon map along the backbone
        # Position 1 to total_carbons
        branch_dict: Dict[int, List[PEBranch]] = {}
        for b in self.branches:
            branch_dict.setdefault(b.position, []).append(b)

        smiles_parts = []
        for pos in range(1, self.total_carbons + 1):
            smiles_parts.append("C")
            if pos in branch_dict:
                for b in branch_dict[pos]:
                    smiles_parts.append(b.to_smiles_substituent())

        return "".join(smiles_parts)

    def to_selfies(self) -> str:
        """Converts the constructed PE SMILES string into SELFIES string."""
        smiles = self.to_smiles()
        if sf is None:
            # Fallback representation when selfies library is not installed
            return f"[C]{smiles}"
        try:
            return sf.encoder(smiles)
        except Exception:
            clean_smiles = smiles.replace("@", "")
            return sf.encoder(clean_smiles)

    @classmethod
    def from_selfies(cls, selfies_str: str) -> "PERepresentation":
        """Reconstructs PE representation from a SELFIES string."""
        if sf is None:
            smiles = selfies_str.replace("[C]", "")
            return cls.from_smiles(smiles)
        try:
            smiles = sf.decoder(selfies_str)
        except Exception:
            smiles = selfies_str.replace("[C]", "")
        return cls.from_smiles(smiles)

    @classmethod
    def from_smiles(cls, smiles: str) -> "PERepresentation":
        """Estimates PE parameters (dp, branches) from a given SMILES string."""
        # Simple heuristic decoder counting carbons and branch points
        total_c = smiles.count("C") + smiles.count("c")
        dp = max(2, total_c // 2)
        branches = []
        
        # Detect branches denoted by parenthesis in SMILES
        import re
        branch_matches = re.findall(r'\(C+?\)', smiles)
        pos_counter = 2
        for match in branch_matches:
            branch_len = match.count("C")
            branches.append(PEBranch(position=pos_counter, length=branch_len))
            pos_counter += 3
            
        return cls(dp=dp, branches=branches)
