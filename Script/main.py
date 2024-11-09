import os
import psutil
import polars as pl
from natsort import natsorted

# Import custom plot functions from separate modules
from Figure1 import *
from Figure2 import *
from Figure3 import *
from Figure4 import *

# Define base directory as the directory containing this script
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "Data")

# Ensure that the data directory exists
if not os.path.exists(DATA_DIR):
    raise FileNotFoundError(f"Data directory not found at {DATA_DIR}")

# Utility functions for system information
def system_info():
    # Number of CPU cores
    num_cores = os.cpu_count()
    print(f"Number of CPU cores: {num_cores}")

    # Amount of RAM
    ram_info = psutil.virtual_memory()
    total_ram = ram_info.total / (1024 ** 3)  # Convert bytes to GB
    available_ram = ram_info.available / (1024 ** 3)  # Convert bytes to GB
    print(f"Total RAM: {total_ram:.2f} GB")
    print(f"Available RAM: {available_ram:.2f} GB")
    print(f'Working Directory: {BASE_DIR}')
    print(f'Data Directory: {DATA_DIR}')

# Plotting functions
def plot_figure1():
    # Figure 1A
    sum_alltype = pl.read_csv(os.path.join(DATA_DIR, "Variant.bed"), separator='\t', has_header=True).with_columns(Count = 1).group_by('col_14','col_8').agg(pl.col('Count').sum())
    sorted_value = natsorted(sum_alltype['col_8'].unique())
    plot_fig1a = Fraction_plot(sum_alltype, sorted_value)

    # Figure 1B & Supplementary S1 A-B
    chr_map = pl.read_csv(os.path.join(DATA_DIR, 'Genome_text.tsv'), separator="\t", has_header=False)
    vcf_10kb_df = pl.read_csv(os.path.join(DATA_DIR, '10Kb_window_Variant_Count.bed'), separator="\t", has_header=False)
    sorted_value = natsorted(Data_processing(vcf_10kb_df, chr_map)['column_3_right'].unique())
    plot_fig1b = Plot_VarChr(Data_processing(vcf_10kb_df, chr_map), sorted_value, save_name='Figure1B')

    # Homozygous Variant Plot (Supplementary Figure S1A)
    homo_df = pl.read_csv(os.path.join(DATA_DIR, "10Kb_window_Variant_Count_SNP_InDel_Homozygous.bed"), separator="\t", has_header=False)
    plot_figS1a = Plot_VarChr(Data_processing(homo_df, chr_map), sorted_value, save_name='Figure_S1A', title='Homozygous variant occurrence per chromosome with resolution of 10Kb')

    # Heterozygous Variant Plot (Supplementary Figure S1B)
    hetero_df = pl.read_csv(os.path.join(DATA_DIR, "10Kb_window_Variant_Count_SNP_InDel_Heterozygous.bed"), separator="\t", has_header=False)
    plot_figS1b = Plot_VarChr(Data_processing(hetero_df, chr_map), sorted_value, save_name='Figure_S1B', title='Heterozygous variant occurrence per chromosome with resolution of 10Kb')

def plot_figure2():
    # Figure 2
    chr_map = pl.read_csv(os.path.join(DATA_DIR, 'Genome_text.tsv'), separator="\t", has_header=False)
    df_500kbp = pl.read_csv(os.path.join(DATA_DIR, "500Kb_window_Variant_Count.bed"), separator="\t", has_header=False)
    centromere_df = pl.read_csv(os.path.join(DATA_DIR, "T2T-MFA8v1.0.centromere.bed"), separator="\t", has_header=False)
    processed_value = Data_processing_fig2(df_500kbp, chr_map, centromere_df)
    sorted_value = natsorted(processed_value['Chromosome'].unique())
    plot_fig2 = Plot_VarPerChr(processed_value, sorted_value)

def plot_figure3():
    # Figure 3 B-C
    plot_fig3bc = plot_sv_chart(os.path.join(DATA_DIR, "SV_type_Median_Length.tsv"))

def plot_figure4():
    # Figure 4
    centromere_df = pl.read_csv(os.path.join(DATA_DIR, "T2T-MFA8v1.0.centromere.bed"), separator="\t", has_header=False)
    SV_df = pl.read_csv(os.path.join(DATA_DIR, "merged_SV_df.tsv"), separator='\t')
    plot_fig4 = Plot_TrendPerChr(Data_processing_fig4(SV_df, centromere_df))

# Main function
if __name__ == "__main__":
    system_info()
    plot_figure1()
    plot_figure2()
    plot_figure3()
    plot_figure4()
