import pandas as pd
import numpy as np
import argparse

def calculate_e_score(summary_stats_files, trait_groups):
    """
    Calculates the Composite E-score based on MTAG summary statistics.
    
    Args:
        summary_stats_files (list): List of file paths to GWAS/MTAG summary stats.
        trait_groups (dict): Mapping of traits to groups (Static vs Dynamic).
    """
    
    # 1. Load Data and Calculate Mean Chi-Square for weights
    weights = {}
    total_weight_denom = 0
    data_frames = {}
    
    print("Calculating weights based on Chi-Square inflation...")
    for filepath in summary_stats_files:
        df = pd.read_csv(filepath, sep='\t')
        # Assuming columns: SNP, P, Z or CHI2
        if 'CHI2' not in df.columns:
            # Convert P to Chi2 if not present: Chi2 = qchisq(1-P, 1) approx
            from scipy.stats import chi2
            df['CHI2'] = chi2.isf(df['P'], 1)
            
        mean_chi2 = df['CHI2'].mean()
        
        # Formula part: max(0, mean_chi2 - 1)
        adjusted_chi2 = max(0, mean_chi2 - 1)
        
        trait_name = filepath.split('/')[-1].replace('.txt', '')
        weights[trait_name] = adjusted_chi2
        total_weight_denom += adjusted_chi2
        
        data_frames[trait_name] = df.set_index('SNP')

    # 2. Compute E-score for each SNP
    print("Computing E-score per SNP...")
    
    # Get common SNPs
    common_snps = set.intersection(*[set(df.index) for df in data_frames.values()])
    
    e_score_results = []
    
    for snp in common_snps:
        e_score_val = 0
        
        for trait_name, df in data_frames.items():
            snp_data = df.loc[snp]
            p_val = snp_data['P']
            
            # Weight calculation: (Sum_j max(0, chi2-1)) / Total_Sum
            # Note: Simplification here assumes calculating per trait contribution
            w = weights[trait_name] / total_weight_denom
            
            # Formula part: - sum ( weight * log10(P) )
            if p_val > 0:
                log_p = np.log10(p_val)
                e_score_val += -1 * (w * log_p)
            
        e_score_results.append({'SNP': snp, 'E_score': e_score_val})
        
    return pd.DataFrame(e_score_results)

if __name__ == "__main__":
    # Example usage boilerplate would go here
    pass
