import os
import gzip
import pandas as pd
from plot import plot_qq, plot_manhattan, plot_pvalue_histogram, plot_effect_size_distribution, plot_volcano, assess_population_stratification
from calculations import ld_score_regression, display_top_snps, calculate_power, calculate_column_statistics, calculate_minor_allele_frequency, calculate_genomic_inflation_factor, calculate_heritability, calculate_effective_sample_size

file_path = 'example_summary_statistics.gz'  # Replace with your .gz file path





def read_gz_file(file_path, n_lines=10):
    required_columns = {'CHR', 'BP', 'ID', 'A1', 'A2', 'BETA', 'SE', 'P', 'INFO', 'NCAS', 'NCON'}
    log_content = []
    
    with gzip.open(file_path, 'rt') as file:
        header = file.readline().strip().split()
        missing_columns = required_columns - set(header)
        
        if missing_columns:
            log_content.append(f"Error: Missing columns: {', '.join(missing_columns)}")
            print("\n".join(log_content))
            return
        
        log_content.append("Header columns found, printing first few lines:")
        log_content.append("\t".join(header))
        
        data = []
        for i, line in enumerate(file):
            if i < n_lines:
                log_content.append(line.strip())
            data.append(line.strip().split())
    
    df = pd.DataFrame(data, columns=header)
    
    numeric_cols = ['CHR', 'BP', 'BETA', 'SE', 'P', 'INFO', 'NCAS', 'NCON']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    
    log_content.append(str(calculate_column_statistics(df, numeric_cols)))
    log_content.append(str(calculate_minor_allele_frequency(df)))
    log_content.append(str(calculate_genomic_inflation_factor(df)))
    log_content.append(str(calculate_heritability(df)))
    log_content.append(str(calculate_effective_sample_size(df)))
    log_content.append(str(ld_score_regression(df)))
    log_content.append(str(display_top_snps(df)))
    log_content.append(str(calculate_power(df)))
    log_content.append(str(assess_population_stratification(df, file_path)))

    plot_qq(df['P'].dropna(), file_path)
    plot_manhattan(df, file_path)
    plot_pvalue_histogram(df['P'].dropna(), file_path)
    plot_effect_size_distribution(df['BETA'].dropna(), file_path)
    plot_volcano(df.dropna(subset=['BETA', 'P']), file_path)
        
    output_dir = os.path.splitext(file_path)[0]
    os.makedirs(output_dir, exist_ok=True)
    
    with open(os.path.join(output_dir, 'log.txt'), 'w', encoding='utf-8') as log_file:
        log_file.write("\n".join(log_content))

if __name__ == "__main__":
    read_gz_file(file_path, n_lines=10)
