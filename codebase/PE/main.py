"""
Main CLI Runner for Polyethylene (PE) Genetic Algorithm Optimization.

Executes the PE structural evolution loop, evaluates thermal degradation index & Mn,
and exports reporting outputs to CSV and RDKit grid images.
"""

import argparse
import sys
import os
from typing import Optional

# Ensure codebase root is in Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from codebase.PE.ga_generator import PEGAGenerator
from codebase.PE.fitness import PEFitnessEvaluator
from codebase.PE.visualizer import PEVisualizer

def run_pe_optimization(
    generations: int = 10,
    pop_size: int = 15,
    target_mn: float = 1000.0,
    output_dir: str = "output_pe"
):
    """Runs the main PE optimization loop."""
    print("=" * 70)
    print(f"Starting Polyethylene (PE) Genetic Algorithm Optimization")
    print(f"Generations: {generations} | Population: {pop_size} | Target Mn: {target_mn} g/mol")
    print("=" * 70)

    evaluator = PEFitnessEvaluator(target_mn=target_mn)
    ga = PEGAGenerator(population_size=pop_size, evaluator=evaluator)
    
    population = ga.initialize_population()
    best_overall_score = -999.0
    best_overall_candidate = None

    for gen in range(1, generations + 1):
        population, scored_pop = ga.step_generation(population)
        top_candidate = scored_pop[0]
        top_props = top_candidate[1]
        
        score = top_props.get("fitness_score", 0.0)
        mn = top_props.get("mn", 0.0)
        risk = top_props.get("degradation_risk_index", 0.0)
        tm = top_props.get("estimated_tm_c", 0.0)
        
        if score > best_overall_score:
            best_overall_score = score
            best_overall_candidate = top_candidate

        print(f"[Gen {gen:02d}/{generations:02d}] Best Score: {score:.4f} | Mn: {mn:.1f} g/mol | Deg Risk: {risk:.2f} | Est Tm: {tm:.1f}°C")

    print("=" * 70)
    print("Optimization Complete!")
    if best_overall_candidate:
        top_props = best_overall_candidate[1]
        print(f"Top Overall PE SMILES: {top_props.get('canonical_smiles')}")
        print(f"Top Score: {top_props.get('fitness_score'):.4f} | Mn: {top_props.get('mn'):.1f} g/mol | Thermal Risk Index: {top_props.get('degradation_risk_index'):.2f}")

    # Generate CSV & Visualization Report
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, "pe_optimization_results.csv")
    img_path = os.path.join(output_dir, "top_pe_candidates.png")
    
    final_scored_pop = ga.evaluate_population(population)
    df = PEVisualizer.create_dataframe(final_scored_pop)
    PEVisualizer.save_dataframe_csv(df, csv_path)
    PEVisualizer.draw_top_molecules(final_scored_pop, img_path, top_k=6)
    
    print(f"\nSaved CSV report to: {csv_path}")
    print(f"Saved RDKit top candidate structures to: {img_path}")
    print("=" * 70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Polyethylene (PE) GA Optimization Runner")
    parser.add_argument("--generations", type=int, default=10, help="Number of GA generations")
    parser.add_argument("--pop_size", type=int, default=15, help="Population size per generation")
    parser.add_argument("--target_mn", type=float, default=1000.0, help="Target Number Average Molecular Weight (g/mol)")
    parser.add_argument("--output_dir", type=str, default="output_pe", help="Output directory for reports")
    
    args = parser.parse_args()
    run_pe_optimization(
        generations=args.generations,
        pop_size=args.pop_size,
        target_mn=args.target_mn,
        output_dir=args.output_dir
    )
