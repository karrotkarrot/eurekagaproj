"""
Unit tests for Polyethylene (PE) optimization package in codebase/PE.
"""

import unittest
import os
import sys

# Ensure codebase root is in Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

from codebase.PE.representation import PERepresentation, PEBranch
from codebase.PE.validator import PEValidator
from codebase.PE.fitness import PEFitnessEvaluator
from codebase.PE.ga_generator import PEGAGenerator
from codebase.PE.visualizer import PEVisualizer

class TestPERepresentation(unittest.TestCase):
    def test_pe_smiles_generation(self):
        pe = PERepresentation(dp=5, branches=[PEBranch(position=3, length=1, bend_direction="straight")])
        smiles = pe.to_smiles()
        self.assertTrue(PEValidator.validate_smiles(smiles))
        self.assertIn("(C)", smiles)

    def test_pe_selfies_conversion(self):
        pe = PERepresentation(dp=4, branches=[PEBranch(position=2, length=2)])
        selfies_str = pe.to_selfies()
        self.assertIsNotNone(selfies_str)
        self.assertTrue(len(selfies_str) > 0)

class TestPEValidator(unittest.TestCase):
    def test_validation_and_canonicalization(self):
        valid_smiles = "CCCCCCCC"
        is_valid, canon = PEValidator.validate_and_canonicalize(valid_smiles)
        self.assertTrue(is_valid)
        self.assertEqual(canon, "CCCCCCCC")

        invalid_smiles = "C(C)(C)(C)(C)C"  # Hypervalent carbon
        is_valid, canon = PEValidator.validate_and_canonicalize(invalid_smiles)
        self.assertFalse(is_valid)

class TestPEFitnessEvaluator(unittest.TestCase):
    def test_fitness_scoring(self):
        evaluator = PEFitnessEvaluator(target_mn=500.0)
        res = evaluator.evaluate_properties("CCCCCCCCCCC")
        self.assertTrue(res["valid"])
        self.assertGreater(res["mn"], 100.0)
        self.assertGreaterEqual(res["fitness_score"], 0.0)

class TestPEGAGenerator(unittest.TestCase):
    def test_ga_step(self):
        ga = PEGAGenerator(population_size=6)
        pop = ga.initialize_population()
        self.assertEqual(len(pop), 6)
        
        next_gen, scored = ga.step_generation(pop)
        self.assertEqual(len(next_gen), 6)
        self.assertTrue(len(scored) > 0)

if __name__ == "__main__":
    unittest.main()
