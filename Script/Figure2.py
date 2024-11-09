import altair as alt
alt.data_transformers.enable("vegafusion")
import polars as pl
## -------- ##

## FIGURE 2 ##
# Use Centromere Length from This Paper (https://doi.org/10.1101/2024.04.07.588379)
# Zhang, S. et al. Comparative genomics of macaques and integrated insights into genetic variation and population history. bioRxiv (2024).
# https://github.com/zhang-shilong/T2T-MFA8

def Data_processing_fig2(df, chr_map, centromere_df):
    
    """
    Function to process data for Figure 2, including joining with chromosome and centromere information and calculating relevant metrics.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data containing genomic variant information.

    chr_map : polars.DataFrame
        A mapping dataframe used to add chromosome information to the main data.

    centromere_df : polars.DataFrame
        A dataframe containing centromere location information, including columns for chromosome, start, and end positions.

    Returns: proc_df
    --------
    polars.DataFrame
        A processed dataframe that includes normalized density, percentage density, and mid-point positions for variant data.

    Notes:
    ------
    - The centromere data is transformed to megabase pairs for easier visualization.
    - The function calculates metrics such as normalized variant density, percentage density, and general positional information.
    - Chromosome names are formatted to maintain a consistent ordering format.

    Example Usage:
    --------------
    >>> processed_df_fig2 = Data_processing_fig2(df, chr_map, centromere_df)
    
    """
    
    centromere_df = centromere_df.rename({'column_1':'Chromosome','column_2':'Start','column_3':'End'})\
    .with_columns(Start = pl.col('Start')/10**6,
                 End = pl.col('End')/10**6)
    
    cleaned_df = df.join(chr_map,'column_1').filter(pl.col('column_3_right')!="MT").with_columns(pl.col('column_2')/10**6, pl.col('column_3')/10**6).select(pl.exclude("column_2_right","column_1")).rename({"column_3_right":"Chromosome"})
    
    joined_df = cleaned_df.join(centromere_df,'Chromosome')
    proc_df = joined_df.with_columns(Norm_den=pl.col('column_4')/500000,
                Percent_den=pl.col('column_4')*100/500000,
                      general_pos=(pl.col('column_2')+pl.col('column_3'))/2)
    return proc_df.sort('Chromosome').with_columns(pl.col('Chromosome').str.replace(r"chr(\d)\b", r"chr0$1"))


def Plot_VarPerChr(df, sorted_value, save_name='Figure2', title="Variant Distribution from Chromosomal End to Centromeric Region at 500 Kbp Resolution", x='general_pos:Q', y='Percent_den:Q'):
    
    """
    Function to create a facet plot of variant distribution from chromosomal ends to centromeric regions.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data containing variant information, with columns for chromosomal positions, centromere start/end, and density metrics.

    sorted_value : list
        A list of chromosome labels to specify the order of chromosomes in the plot.

    save_name : str, optional (default: 'Figure2')
        The filename used to save the generated plot. The plot will be saved in both .svg and .png formats.

    title : str, optional (default: 'Variant Distribution from Chromosomal End to Centromeric Region at 500 Kbp Resolution')
        The title of the plot.

    x : str, optional (default: 'general_pos:Q')
        The column to be used for the x-axis, representing the position along the chromosome in megabase pairs.

    y : str, optional (default: 'Percent_den:Q')
        The column to be used for the y-axis, representing the percentage density of variants.

    Returns:
    --------
    altair.Chart
        The generated Altair facet plot, displaying variant distribution across chromosomes.

    Notes: combined_plot
    ------
    - The plot uses a combination of scatter and rectangle marks to visualize variant distribution and centromere regions.
    - The generated plot uses category20b colors and is saved at a resolution suitable for publication.

    Example Usage:
    --------------
    >>> plot_fig2 = Plot_VarPerChr(df, sorted_value=['1', '2', '3', 'X', 'Y'], save_name='VariantDistributionPlot')
    
    """
    
    color_scale = alt.Scale(scheme="category20b")

    plot = alt.Chart(df,
            title=title
            ).encode(
        x = alt.X(x,title='Position (Mbp)'),
        y = alt.Y(y,title='Percentage of Variant'),
        color=alt.Color('Chromosome',scale=color_scale,sort=sorted_value,title=None)
    ).properties(
        width=200,
        height=100
    )

    combined_plot = alt.layer(
        plot.mark_circle(size=40),
        plot.mark_rect(color='', fill='', stroke='grey', strokeWidth=1.4, strokeDash=[2, 2]).encode(
            x='Start:Q',
            x2='End:Q',
            y=alt.value(0),  # 0 pixels from top
        )
    ).facet(
        "Chromosome:N",
        columns=3
    ).configure_header(
        titleFontSize=16,  # Change this value to adjust the font size of the facet titles
        labelFontSize=14   # Optional: Adjusts the label font size of each facet as well
    ).resolve_axis(
        x='independent',
        y='independent'
    ).resolve_scale(
        x='independent', 
        y='independent'
    )
    
    # Save
    combined_plot.save(f'{save_name}.svg', engine='vl-convert', ppi=300)
    combined_plot.save(f'{save_name}.png', engine='vl-convert', ppi=300)
    
    return combined_plot




