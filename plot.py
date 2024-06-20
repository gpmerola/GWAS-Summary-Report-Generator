import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def plot_qq(p_values, file_path):
    output_dir = os.path.splitext(file_path)[0]
    sorted_pvals = np.sort(p_values)
    expected = -np.log10(np.arange(1, len(sorted_pvals) + 1) / (len(sorted_pvals) + 1))
    observed = -np.log10(sorted_pvals)
    
    plt.figure(figsize=(6, 6))
    plt.plot(expected, observed, 'o', markersize=2)
    plt.plot([0, max(expected)], [0, max(expected)], 'r--')
    plt.xlabel('Expected -log10(p)')
    plt.ylabel('Observed -log10(p)')
    plt.title('QQ Plot')
    plt.savefig(os.path.join(output_dir, 'qq_plot.png'))
    plt.close()

def plot_manhattan(df, file_path):
    output_dir = os.path.splitext(file_path)[0]
    df['-log10(P)'] = -np.log10(df['P'])
    df = df.dropna(subset=['CHR', 'BP', '-log10(P)'])
    df = df[np.isfinite(df['-log10(P)'])]
    df['CHR'] = df['CHR'].astype(int)
    chr_dict = {chrom: index for index, chrom in enumerate(sorted(df['CHR'].unique()))}
    df['chr_index'] = df['CHR'].map(chr_dict)
    df['normalized_bp'] = df.groupby('CHR')['BP'].transform(lambda x: (x - x.min()) / (x.max() - x.min()))

    fig, ax = plt.subplots(figsize=(18, 6))
    colors = ['b', 'r']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df.groupby('CHR')):
        ax.scatter(group['chr_index'] + group['normalized_bp'], group['-log10(P)'], color=colors[num % len(colors)], s=10)
        x_labels.append(name)
        x_labels_pos.append(chr_dict[name] + 0.5)
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, fontsize=8)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Manhattan Plot')
    plt.savefig(os.path.join(output_dir, 'manhattan_plot.png'))
    plt.close()

def plot_pvalue_histogram(p_values, file_path):
    output_dir = os.path.splitext(file_path)[0]
    plt.figure(figsize=(8, 6))
    plt.hist(p_values, bins=50, edgecolor='k', alpha=0.7)
    plt.xlabel('P-value')
    plt.ylabel('Frequency')
    plt.title('P-value Distribution')
    plt.savefig(os.path.join(output_dir, 'pvalue_histogram.png'))
    plt.close()

def plot_effect_size_distribution(effect_sizes, file_path):
    output_dir = os.path.splitext(file_path)[0]
    plt.figure(figsize=(8, 6))
    plt.hist(effect_sizes, bins=50, edgecolor='k', alpha=0.7)
    plt.xlabel('Effect Size (BETA)')
    plt.ylabel('Frequency')
    plt.title('Effect Size Distribution')
    plt.savefig(os.path.join(output_dir, 'effect_size_distribution.png'))
    plt.close()

def plot_volcano(df, file_path):
    output_dir = os.path.splitext(file_path)[0]
    df['-log10(P)'] = -np.log10(df['P'])
    
    plt.figure(figsize=(10, 8))
    plt.scatter(df['BETA'], df['-log10(P)'], alpha=0.5)
    plt.xlabel('Effect Size (BETA)')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot')
    plt.savefig(os.path.join(output_dir, 'volcano_plot.png'))
    plt.close()

def assess_population_stratification(df, file_path, n_components=2):
    output_dir = os.path.splitext(file_path)[0]
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(df[['NCAS', 'NCON']].dropna())
    
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.5)
    plt.xlabel('PCA1')
    plt.ylabel('PCA2')
    plt.title('Population Stratification Assessment')
    plt.savefig(os.path.join(output_dir, 'population_stratification.png'))
    plt.close()
