#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: calc_E_score.py
Description:
    Calculates the composite E-score by integrating summary statistics from 
    multiple traits (e.g., Static weights + Dynamic growth parameters).
    
    The E-score is a weighted sum of -log10(P-values), where weights are 
    determined by the mean Chi-Square statistic of each trait, emphasizing 
    traits with stronger genomic inflation (higher heritability/signal).

Methodology:
    E(i) = Sum_k [ w_k * (-log10(P_{k,i})) ]
    where w_k = max(0, mean(Chi2_k) - 1) / Sum_all(max(0, mean(Chi2) - 1))

Usage:
    python calc_E_score.py --input data/gwas_results/*.assoc.txt --output output/escore_results.txt

Dependencies:
    pandas, numpy, scipy
"""

import argparse
import sys
import pandas as pd
import numpy as np
from scipy.stats import chi2
from functools import reduce
from pathlib import Path

# =============================================================================
# Utility Functions
# =============================================================================

def parse_arguments():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate Composite E-score from GWAS Summary Statistics."
    )
    parser.add_argument(
        '--input', 
        nargs='+', 
        required=True, 
        help='List of input summary statistic files (e.g., .assoc.txt or .summary).'
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help='Path to save the output E-score file.'
    )
    parser.add_argument(
        '--p_col', 
        default='P', 
        help='Column name for P-value in input files (default: "P").'
    )
    parser.add_argument(
        '--snp_col', 
        default='SNP', 
        help='Column name for SNP ID in input files (default: "SNP").'
    )
    return parser.parse_args()

def load_and_calculate_weight(filepath, p_col_name='P'):
    """
    Loads a single summary stat file and calculates its weight contribution.
    
    Returns:
        tuple: (trait_name, DataFrame, weight_numerator)
    """
    try:
        # Read file (assuming tab-separated, auto-detect header)
        df = pd.read_csv(filepath, sep=None, engine='python')
        
        # Check required columns
        if p_col_name not in df.columns:
            raise ValueError(f"Column '{p_col_name}' not found in {filepath}")
            
        # Calculate Chi-Square statistics if not present
        if 'CHI2' not in df.columns:
            # Convert P to Chi2: ISF is Inverse Survival Function (Inverse CDF)
            # chi2.isf(P, df=1) approximates the Z^2 statistic
            df['CHI2'] = chi2.isf(df[p_col_name], 1)

        # Calculate Mean Chi-Square Inflation
        mean_chi2 = df['CHI2'].mean()
        
        # Determine Weight Numerator: max(0, mean_chi2 - 1)
        # We subtract 1 because the expected mean Chi2 under null hypothesis is 1.
        weight_numerator = max(0, mean_chi2 - 1)
        
        trait_name = Path(filepath).stem  # Extract filename without extension
        print(f"  [Loaded] {trait_name}: Mean Chi2 = {mean_chi2:.4f}, Weight Score = {weight_numerator:.4f}")
        
        return trait_name, df, weight_numerator

    except Exception as e:
        print(f"Error processing file {filepath}: {e}", file=sys.stderr)
        sys.exit(1)

# =============================================================================
# Main Logic
# =============================================================================

def main():
    args = parse_arguments()
    
    print("-" * 60)
    print(f"Starting E-score calculation for {len(args.input)} traits...")
    print("-" * 60)

    # 1. Load Data and Calculate Weights
    # -----------------------------------------------------------
    loaded_data = []
    total_weight_denom = 0.0
    
    for filepath in args.input:
        trait_name, df, w_num = load_and_calculate_weight(filepath, args.p_col)
        
        # Filter strictly necessary columns to save memory
        df = df[[args.snp_col, args.p_col]].copy()
        
        # Rename P-value column to avoid collision after merge
        df.rename(columns={args.p_col: f"P_{trait_name}"}, inplace=True)
        
        loaded_data.append({
            'name': trait_name,
            'df': df,
            'weight_score': w_num
        })
        total_weight_denom += w_num

    if total_weight_denom == 0:
        print("Error: Total weight denominator is 0. No inflation detected in any trait.", file=sys.stderr)
        sys.exit(1)

    # 2. Vectorized E-score Calculation
    # -----------------------------------------------------------
    print("\nMerging datasets to identify common SNPs...")
    
    # Use functools.reduce to merge all DataFrames on SNP column (Inner Join)
    # This ensures we only calculate E-scores for SNPs present in ALL summary stats
    merged_df = reduce(
        lambda left, right: pd.merge(left['df'], right['df'], on=args.snp_col), 
        loaded_data
    )
    
    print(f"  Common SNPs identified: {len(merged_df):,}")
    print("Computing composite E-scores (Vectorized)...")

    # Initialize E-score column
    merged_df['E_score'] = 0.0
    
    # Vectorized accumulation: Sum( -log10(P) * Normalized_Weight )
    for item in loaded_data:
        trait_name = item['name']
        w_score = item['weight_score']
        
        # Calculate Normalized Weight
        normalized_w = w_score / total_weight_denom
        
        p_col = f"P_{trait_name}"
        
        # Add contribution to E-score
        # E += -log10(P) * w
        # Use np.maximum to avoid log(0) error, though P-values should be > 0
        p_values = np.maximum(merged_df[p_col], 1e-300) 
        merged_df['E_score'] += -np.log10(p_values) * normalized_w

    # 3. Save Results
    # -----------------------------------------------------------
    # Select only SNP and E_score for output (clean format)
    output_df = merged_df[[args.snp_col, 'E_score']].sort_values(by='E_score', ascending=False)
    
    print(f"\nSaving results to: {args.output}")
    output_df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    
    print("Done.")

if __name__ == "__main__":
    main()
