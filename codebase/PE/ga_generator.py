"""
PE Genetic Algorithm (GA) Generator Module.

Uses SELFIES representation for robust genetic operations (crossover & mutation)
to optimize PE polymer structures against thermal degradation and maximize Mn.
"""

import random
from typing import List, Dict, Any, Tuple
try:
    import selfies as sf
except ImportError:
    sf = None
from codebase.PE.representation import PERepresentation, PEBranch
from codebase.PE.validator import PEValidator
from codebase.PE.fitness import PEFitnessEvaluator

class PEGAGenerator:
    """Genetic Algorithm Optimizer for Polyethylene (PE) Candidate Generation."""

    def __init__(
        self,
        population_size: int = 20,
        mutation_rate: float = 0.3,
        crossover_rate: float = 0.7,
        evaluator: Optional[PEFitnessEvaluator] = None
    ):
        """
        :param population_size: Number of PE individuals per generation.
        :param mutation_rate: Probability of mutation per individual.
        :param crossover_rate: Probability of crossover between parent pairs.
        :param evaluator: PEFitnessEvaluator instance.
        """
        self.pop_size = max(4, population_size)
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.evaluator = evaluator if evaluator is not None else PEFitnessEvaluator()
        self.evaluated_cache: Dict[str, Dict[str, Any]] = {}

    def initialize_population(self) -> List[PERepresentation]:
        """Seeds an initial population of diverse PE structures (linear & branched)."""
        population = []
        for i in range(self.pop_size):
            dp = random.randint(5, 25)
            num_branches = random.randint(0, 4)
            branches = []
            for _ in range(num_branches):
                pos = random.randint(2, dp * 2 - 1)
                length = random.choice([1, 2, 4])  # Methyl, Ethyl, or Butyl branches
                bend = random.choice(["left", "right", "straight"])
                branches.append(PEBranch(pos, length, bend))
            
            pe = PERepresentation(dp=dp, branches=branches)
            population.append(pe)
        return population

    def evaluate_population(self, population: List[PERepresentation]) -> List[Tuple[PERepresentation, Dict[str, Any]]]:
        """Evaluates fitness for each PE individual in the population using canonical caching."""
        scored_pop = []
        for pe in population:
            smiles = pe.to_smiles()
            canon = PEValidator.canonicalize(smiles)
            if canon in self.evaluated_cache:
                props = self.evaluated_cache[canon]
            else:
                props = self.evaluator.evaluate_properties(smiles)
                if canon:
                    self.evaluated_cache[canon] = props
            scored_pop.append((pe, props))
        
        # Sort population by fitness score descending
        scored_pop.sort(key=lambda x: x[1]["fitness_score"], reverse=True)
        return scored_pop

    def mutate_selfies(self, selfies_str: str) -> str:
        """Mutates a SELFIES string by modifying, inserting, or deleting SELFIES symbols."""
        if sf is None:
            # Fallback simple string mutation when selfies library is not installed
            if "[C]" in selfies_str:
                return selfies_str.replace("[C]", "[C][C]", 1)
            return selfies_str + "[C]"

        tokens = list(sf.split_selfies(selfies_str))
        if not tokens:
            return selfies_str
        
        mutation_type = random.choice(["modify", "insert", "delete"])
        if mutation_type == "modify" and len(tokens) > 0:
            idx = random.randint(0, len(tokens) - 1)
            tokens[idx] = random.choice(["[C]", "[Branch1]", "[Ring1]", "[=C]"])
        elif mutation_type == "insert":
            idx = random.randint(0, len(tokens))
            tokens.insert(idx, "[C]")
        elif mutation_type == "delete" and len(tokens) > 1:
            idx = random.randint(0, len(tokens) - 1)
            tokens.pop(idx)
        
        mutated = "".join(tokens)
        try:
            smiles = sf.decoder(mutated)
            if PEValidator.validate_smiles(smiles):
                return mutated
        except Exception:
            pass
        return selfies_str

    def crossover_selfies(self, parent1_selfies: str, parent2_selfies: str) -> Tuple[str, str]:
        """Performs single-point crossover on two parent SELFIES strings."""
        if sf is None:
            mid1 = len(parent1_selfies) // 2
            mid2 = len(parent2_selfies) // 2
            return (parent1_selfies[:mid1] + parent2_selfies[mid2:], parent2_selfies[:mid2] + parent1_selfies[mid1:])

        t1 = list(sf.split_selfies(parent1_selfies))
        t2 = list(sf.split_selfies(parent2_selfies))
        
        if len(t1) < 2 or len(t2) < 2:
            return parent1_selfies, parent2_selfies
        
        cut1 = random.randint(1, len(t1) - 1)
        cut2 = random.randint(1, len(t2) - 1)
        
        child1_tokens = t1[:cut1] + t2[cut2:]
        child2_tokens = t2[:cut2] + t1[cut1:]
        
        c1_str = "".join(child1_tokens)
        c2_str = "".join(child2_tokens)
        
        return c1_str, c2_str

    def step_generation(
        self, population: List[PERepresentation]
    ) -> Tuple[List[PERepresentation], List[Tuple[PERepresentation, Dict[str, Any]]]]:
        """Advances the population by 1 generation (Evaluate -> Select -> Crossover -> Mutate)."""
        scored_pop = self.evaluate_population(population)
        
        # Elitism: Retain top 20%
        num_elites = max(1, int(self.pop_size * 0.2))
        next_gen: List[PERepresentation] = [item[0] for item in scored_pop[:num_elites]]
        
        # Selection pool (Tournament selection)
        def select_parent():
            candidates = random.sample(scored_pop, k=min(3, len(scored_pop)))
            candidates.sort(key=lambda x: x[1]["fitness_score"], reverse=True)
            return candidates[0][0]

        while len(next_gen) < self.pop_size:
            p1 = select_parent()
            p2 = select_parent()
            
            p1_selfies = p1.to_selfies()
            p2_selfies = p2.to_selfies()
            
            # Crossover
            if random.random() < self.crossover_rate:
                c1_selfies, c2_selfies = self.crossover_selfies(p1_selfies, p2_selfies)
            else:
                c1_selfies, c2_selfies = p1_selfies, p2_selfies
            
            # Mutation
            if random.random() < self.mutation_rate:
                c1_selfies = self.mutate_selfies(c1_selfies)
            if random.random() < self.mutation_rate:
                c2_selfies = self.mutate_selfies(c2_selfies)
                
            # Convert back to PERepresentation
            child1 = PERepresentation.from_selfies(c1_selfies)
            child2 = PERepresentation.from_selfies(c2_selfies)
            
            next_gen.append(child1)
            if len(next_gen) < self.pop_size:
                next_gen.append(child2)
                
        return next_gen, scored_pop
