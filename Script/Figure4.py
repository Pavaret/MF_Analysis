import altair as alt
alt.data_transformers.enable("vegafusion")
import polars as pl
## -------- ##
## FIGURE 4 ##
## Use Centromere Length from Newly Published Paper (https://doi.org/10.1101/2024.04.07.588379)
# Zhang, S. et al. Comparative genomics of macaques and integrated insights into genetic variation and population history. bioRxiv (2024).

# https://github.com/zhang-shilong/T2T-MFA8



def Data_processing_fig4(df, centromere_df):
    """
    Function to process data for Figure 4, including joining with centromere information and calculating relevant metrics.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data containing genomic variant information.

    centromere_df : polars.DataFrame
        A dataframe containing centromere location information, including columns for chromosome, start, and end positions.

    Returns:
    --------
    polars.DataFrame
        A processed dataframe that includes normalized density, percentage density, and mid-point positions for variant data.

    Notes:
    ------
    - The centromere data is transformed to megabase pairs for easier visualization.
    - The function calculates metrics such as normalized variant density, percentage density, and general positional information.
    - Chromosome names are formatted to maintain a consistent ordering format.
    - The `Type` column is modified to replace occurrences of "BND" with "TRN".

    Example Usage:
    --------------
    >>> processed_df_fig4 = Data_processing_fig4(df, centromere_df)
    
    """
    centromere_df = centromere_df.rename({'column_1':'Chromosome','column_2':'Start','column_3':'End'})\
    .with_columns(Start = pl.col('Start')/10**6,
                 End = pl.col('End')/10**6)
    
    df = df.select(pl.exclude("column_2_right","column_1")).join(centromere_df,'Chromosome').with_columns(Norm_den=pl.col('column_4')/500000, Percent_den=pl.col('column_4')*100/500000, general_pos=(pl.col('column_2')+pl.col('column_3'))/2)
    
    proc_df = df.with_columns(pl.col('Type').str.replace("BND", "TRN"))
    
    return proc_df.sort('Chromosome').with_columns(pl.col('Chromosome').str.replace(r"chr(\d)\b", r"chr0$1"))


def Plot_TrendPerChr(df, save_name='Figure4', title="Variant Distribution from Chromosomal End to Centromeric Region at 500 Kbp Resolution", x='column_2:Q', y='column_4:Q'):
    """
    Function to create a facet plot of variant trends per chromosome, showing variant distribution from chromosomal ends to centromeric regions.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data containing variant information, with columns for chromosomal positions, centromere start/end, and density metrics.

    save_name : str, optional (default: 'Figure4')
        The filename used to save the generated plot. The plot will be saved in both .svg and .png formats.

    title : str, optional (default: 'Variant Distribution from Chromosomal End to Centromeric Region at 500 Kbp Resolution')
        The title of the plot.

    x : str, optional (default: 'column_2:Q')
        The column to be used for the x-axis, representing the position along the chromosome in megabase pairs.

    y : str, optional (default: 'column_4:Q')
        The column to be used for the y-axis, representing the percentage density of variants.

    Returns:
    --------
    altair.Chart
        The generated Altair facet plot, displaying variant distribution and trends across chromosomes.

    Notes:
    ------
    - The plot uses a combination of scatter and rectangle marks to visualize variant distribution and centromere regions.
    - A LOESS smoothed line is included to indicate general trends for each variant type.
    - The generated plot uses set1 colors and is saved at a high resolution suitable for publication.

    Example Usage:
    --------------
    >>> plot_fig4 = Plot_TrendPerChr(df, save_name='VariantTrendPlot')
    
    """
    color_scale = alt.Scale(scheme="set1")
    plot = alt.Chart(df,
            title=title
            ).encode(
        x = alt.X(x,title='Position (Mbp)'),
        y = alt.Y(y,title='Percentage of Variant'),
        color=alt.Color('Type:N',
                        scale=color_scale,
                        sort=df['Type'].unique().sort().to_list(),
                        title=None)
    ).properties(
        width=300,
        height=200
    )

    combined_plot = alt.layer(
        plot.mark_circle(size=40,opacity=0.3),
        plot.mark_rect(color='', fill='',stroke='grey', strokeWidth=1.4, strokeDash=[2, 2]).encode(x = 'Start:Q', x2 = 'End:Q', 
                                        y=alt.value(0),  # 0 pixels from top
                                        ),
        plot.transform_loess('column_2', 'column_4', groupby=['Type']).mark_line(size=4.5).encode(opacity = alt.value(0.75))
    ).facet("Chromosome:N",
            columns = 3
    ).configure_title(
        fontSize=20,
        anchor='middle',
        align='left'
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

    # save plot
    combined_plot.save(f'{save_name}.svg', engine="vl-convert", ppi=450)
    combined_plot.save(f'{save_name}.png', engine="vl-convert", ppi=450)

    return combined_plot
