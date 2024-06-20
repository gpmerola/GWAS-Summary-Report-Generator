# GWAS-Summary-Report-Generator

## Overview
This project is designed to read and analyze GWAS (Genome-Wide Association Study) summary statistics from a .gz file. It performs various statistical calculations and generates several plots to visualize the data.

## Structure
main.py: Main script to read the .gz file, perform calculations, and generate plots.
plot.py: Contains functions to generate QQ plot, Manhattan plot, p-value histogram, effect size distribution plot, volcano plot, and population stratification assessment.
calculations.py: Contains functions for LD score regression, displaying top SNPs, calculating power, column statistics, minor allele frequency, genomic inflation factor, heritability, and effective sample size.

## Usage
Ensure you have Python 3 installed.
Install the necessary dependencies (see below).
Place your GWAS summary statistics .gz file in the project directory.
Modify the file_path variable in main.py to point to your .gz file.
Run the main script:
sh
Copia codice
python main.py
## Dependencies
Install the required Python packages using pip. You can create a virtual environment to manage dependencies.

### Dependencies Installation
Create and activate a virtual environment (optional but recommended):

sh
Copia codice
python -m venv venv
source venv/bin/activate   # On Windows use `venv\Scripts\activate`
Install the required packages:

sh
Copia codice
pip install -r requirements.txt
requirements.txt
Copia codice
numpy
pandas
matplotlib
scipy
statsmodels
scikit-learn

## Output
The script will generate the following output files in a directory named after the input file (without the .gz extension):

log.txt: Contains logs and summaries of the calculations.
qq_plot.png: QQ plot of p-values.
manhattan_plot.png: Manhattan plot of p-values.
pvalue_histogram.png: Histogram of p-values.
effect_size_distribution.png: Histogram of effect sizes.
volcano_plot.png: Volcano plot of effect sizes vs. p-values.
population_stratification.png: PCA plot for population stratification assessment.

## Functions
main.py
read_gz_file(file_path, n_lines=10): Reads the .gz file, performs various calculations, generates plots, and logs the results.
plot.py
plot_qq(p_values, file_path): Generates a QQ plot of p-values.
plot_manhattan(df, file_path): Generates a Manhattan plot of p-values.
plot_pvalue_histogram(p_values, file_path): Generates a histogram of p-values.
plot_effect_size_distribution(effect_sizes, file_path): Generates a histogram of effect sizes.
plot_volcano(df, file_path): Generates a volcano plot of effect sizes vs. p-values.
assess_population_stratification(df, file_path, n_components=2): Assesses population stratification using PCA.
calculations.py
ld_score_regression(df): Performs LD score regression.
display_top_snps(df, top_n=10): Displays the top SNPs based on p-values.
calculate_power(df, alpha=0.05, effect_size=0.2): Calculates the statistical power of the study.
calculate_column_statistics(df, numeric_cols): Calculates mean and standard deviation of numeric columns.
calculate_minor_allele_frequency(df): Calculates the minor allele frequency (MAF).
calculate_genomic_inflation_factor(df): Calculates the genomic inflation factor (λ GC).
calculate_heritability(df): Calculates narrow-sense heritability (h²).
calculate_effective_sample_size(df): Calculates the effective sample size (NEFF).
