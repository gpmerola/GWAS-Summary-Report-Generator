import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.power import NormalIndPower

def ld_score_regression(df):
    df['chi2'] = stats.chi2.isf(df['P'], 1)
    ldsc_model = sm.OLS(df['chi2'], sm.add_constant(df['INFO']))
    results = ldsc_model.fit()
    
    output = f"\nLD Score Regression Results:\n{results.summary()}\n"
    ld_h2 = results.params['INFO']
    output += f"\nEstimated Heritability from LD Score Regression: {ld_h2}"
    
    print(output)
    return output

def display_top_snps(df, top_n=10):
    top_snps = df.nsmallest(top_n, 'P')
    output = "\nTop SNPs based on p-values:\n"
    output += top_snps[['CHR', 'BP', 'ID', 'BETA', 'SE', 'P']].to_string(index=False)
    
    print(output)
    return output

def calculate_power(df, alpha=0.05, effect_size=0.2):
    analysis = NormalIndPower()
    power = analysis.solve_power(effect_size=effect_size, nobs1=len(df), alpha=alpha)
    output = f"\nStatistical Power of the Study (alpha={alpha}, effect size={effect_size}): {power}"
    
    print(output)
    return output

def calculate_column_statistics(df, numeric_cols):
    output = "\nColumn Statistics:\n"
    for col in numeric_cols:
        if col in df.columns and col not in ['BP', 'CHR']:
            mean = df[col].mean()
            std = df[col].std()
            output += f"{col} - Mean: {mean}, SD: {std}\n"
    
    print(output)
    return output

def calculate_minor_allele_frequency(df):
    df['MAF'] = (df['NCAS'] + df['NCON']) / (2 * (df['NCAS'] + df['NCON'] + df['NCAS'] + df['NCON']))
    output = f"\nMAF - Mean: {df['MAF'].mean()}, SD: {df['MAF'].std()}"
    
    print(output)
    return output

def calculate_genomic_inflation_factor(df):
    chi_squared = stats.chi2.isf(df['P'], 1)
    lambda_gc = np.median(chi_squared) / stats.chi2.ppf(0.5, 1)
    output = f"\nGenomic Inflation Factor (λ GC): {lambda_gc}"
    
    print(output)
    return output

def calculate_heritability(df):
    df['h2'] = 2 * df['MAF'] * (1 - df['MAF']) * (df['BETA'] ** 2)
    heritability = df['h2'].sum()
    output = f"\nNarrow-sense Heritability (h²): {heritability}"
    
    print(output)
    return output

def calculate_effective_sample_size(df):
    df['NEFF'] = (4 * df['NCAS'] * df['NCON']) / (df['NCAS'] + df['NCON'])
    neff_mean = df['NEFF'].mean()
    output = f"\nEffective Sample Size (NEFF): {neff_mean}"
    
    print(output)
    return output
