# GWAS-Summary-Report-Generator
This project is designed to read and analyze GWAS (Genome-Wide Association Study) summary statistics from a .gz file. It performs various statistical calculations and generates several plots to visualize the data.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Example](#example)
- [License](#license)

## Structure
- `main.py`: Main script to read the `.gz` file, perform calculations, and generate plots.
- `plot.py`: Contains functions to generate QQ plot, Manhattan plot, p-value histogram, effect size distribution plot, volcano plot, and population stratification assessment.
- `calculations.py`: Contains functions for LD score regression, displaying top SNPs, calculating power, column statistics, minor allele frequency, genomic inflation factor, heritability, and effective sample size.

## Usage
1. Ensure you have Python 3 installed.
2. Install the necessary dependencies (see below).
3. Place your GWAS summary statistics `.gz` file in the project directory.
4. Modify the `file_path` variable in `main.py` to point to your `.gz` file.
5. Run the main script:

   ```sh
   python main.py
   ```

The `.gz` file must contain the following columns for the script to work correctly:

- `CHR`: Chromosome number
- `BP`: Base-pair position
- `ID`: SNP identifier
- `A1`: First allele
- `A2`: Second allele
- `BETA`: Effect size estimate
- `SE`: Standard error of the effect size estimate
- `P`: P-value for the association
- `INFO`: Imputation information score
- `NCAS`: Number of cases
- `NCON`: Number of controls

If any of these columns are missing, the script will log an error and stop execution.

   
## Dependencies
Create and activate a virtual environment (optional but recommended):

   ```sh
   python -m venv venv
   source venv/bin/activate   # On Windows use `venv\Scripts\activate`
   ```

Install the required packages:

   ```sh
   pip install -r requirements.txt
   requirements.txt
   Copia codice
   numpy
   pandas
   matplotlib
   scipy
   statsmodels
   scikit-learn
   ```

## Output
The script will generate the following output files in a directory named after the input file (without the .gz extension):

- log.txt: Contains logs and summaries of the calculations.
- qq_plot.png: QQ plot of p-values.
- manhattan_plot.png: Manhattan plot of p-values.
- pvalue_histogram.png: Histogram of p-values.
- effect_size_distribution.png: Histogram of effect sizes.
- volcano_plot.png: Volcano plot of effect sizes vs. p-values.
- population_stratification.png: PCA plot for population stratification assessment.

## Functions
1. main.py
   - read_gz_file(file_path, n_lines=10): Reads the .gz file, performs various calculations, generates plots, and logs the results.
2. plot.py
   - plot_qq(p_values, file_path): Generates a QQ plot of p-values.
   - plot_manhattan(df, file_path): Generates a Manhattan plot of p-values.
   - plot_pvalue_histogram(p_values, file_path): Generates a histogram of p-values.
   - plot_effect_size_distribution(effect_sizes, file_path): Generates a histogram of effect sizes.
   - plot_volcano(df, file_path): Generates a volcano plot of effect sizes vs. p-values.
   - assess_population_stratification(df, file_path, n_components=2): Assesses population stratification using PCA.
3. calculations.py
   - ld_score_regression(df): Performs LD score regression.
   - display_top_snps(df, top_n=10): Displays the top SNPs based on p-values.
   - calculate_power(df, alpha=0.05, effect_size=0.2): Calculates the statistical power of the study.
   - calculate_column_statistics(df, numeric_cols): Calculates mean and standard deviation of numeric columns.
   - calculate_minor_allele_frequency(df): Calculates the minor allele frequency (MAF).
   - calculate_genomic_inflation_factor(df): Calculates the genomic inflation factor (λ GC).
   - calculate_heritability(df): Calculates narrow-sense heritability (h²).
   - calculate_effective_sample_size(df): Calculates the effective sample size (NEFF).

## License
This project is licensed under the MIT License. See the LICENSE file for details.
