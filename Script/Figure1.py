import altair as alt
alt.data_transformers.enable("vegafusion")
import polars as pl
## -------- ##
## FIGURE 1A ##
def Fraction_plot(df, sorted_value, save_name='Figure1A',
                  y='col_8', x='Count', color='col_14', label="{'Exon':'Exon', 'Lnc_RNA':'lnc RNA', 'Pseudogene':'Pseudogene', 'Intron':'Intron', 'Intergene':'Intergene'}[datum.label]"):
    """
    Function to create a bar plot representing the fraction of genetic variants across different genomic regions.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data to be visualized. It should contain columns representing genomic features and counts.

    sorted_value : list
        A list of categories to specify the order in which the y-axis should be sorted.

    save_name : str, optional (default: 'Figure1A')
        The filename used to save the generated plot. The plot will be saved in both .svg and .png formats.

    y : str, optional (default: 'col_8')
        The column name in the dataframe to be used for the y-axis, representing different genomic regions.

    x : str, optional (default: 'Count')
        The column name in the dataframe to be used for the x-axis, representing the count of genetic variants.

    color : str, optional (default: 'col_14')
        The column name in the dataframe to be used for color encoding, representing different genomic regions.

    label : str, optional (default: "{'Exon':'Exon', 'Lnc_RNA':'lnc RNA', 'Pseudogene':'Pseudogene', 'Intron':'Intron', 'Intergene':'Intergene'}[datum.label]")
        A dictionary-style label mapping for translating region labels in the legend.

    Returns: plot
    --------
    None
        The function saves the plot as both .svg and .png files, but does not return any objects.

    Notes:
    ------
    - The plot displays the fraction of genetic variants normalized across different genomic regions.
    - The generated plot uses pastel colors and is saved at a resolution suitable for publication.

    Example Usage:
    --------------
    >>> Fraction_plot(df, sorted_value=['Exon', 'Lnc_RNA', 'Pseudogene', 'Intron', 'Intergene'], save_name='GenomicFractionPlot')
    
    """
    color_scale = alt.Scale(scheme="pastel1")
    
    plot = alt.Chart(df, title="Fraction of Genetic Variants across Different Genomic Regions").mark_bar().encode(
        y=alt.Y(y, title="", sort=sorted_value, axis=alt.Axis(grid=True,ticks=False)),
        x=alt.X(x, title="Fraction [%]").stack("normalize"),
        color=alt.Color(color, scale=color_scale, title="Region",
                        legend=alt.Legend(
                            title="Region",
                            labelExpr=label
                        )
                    )
    ).properties(
        width=500,
        height=300,
    )

    # Save
    plot.save(f'{save_name}.svg', engine='vl-convert', ppi=300)
    plot.save(f'{save_name}.png', engine='vl-convert', ppi=300)
    
    return plot


## FIGURE 1B + SUPPLEMENTARY_FIGURE S1 A-B ##
## DATA PROCESSING 
def Data_processing(df, chr_map):
    """
    Function to process genetic variant data and add normalized variant density.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data containing variant information. It should have columns for chromosome, start, end positions, and other relevant features.

    chr_map : polars.DataFrame
        A mapping dataframe used to join with the main data to add chromosome information.

    Returns:
    --------
    polars.DataFrame
        A processed dataframe with additional columns for normalized variant density and percentage density.

    Notes:
    ------
    - The function computes normalized variant density by dividing the variant count by the genomic length.
    - The final dataframe is filtered to exclude mitochondrial chromosome data.

    Example Usage:
    --------------
    >>> processed_df = Data_processing(df, chr_map)
    
    """
    # Add normalized variant density
    proc_df = df.with_columns(Norm_Var_Den=pl.col('column_4') / (pl.col('column_3') - pl.col('column_2')))

    # Process and filter data
    plot_df = proc_df.with_columns(
        Percent_den=pl.col('Norm_Var_Den') * 100
    ).join(chr_map, 'column_1').filter(pl.col('column_3_right') != "MT")
    
    return plot_df


def Plot_VarChr(df, sorted_value, save_name='Figure1B', y="Percent_den:Q", x='column_3_right', color='column_3_right:N', title='Variant occurrence per chromosome with resolution of 10Kb'):
    """
    Function to create a scatter plot representing the percentage of variant occurrences per chromosome.

    Parameters:
    -----------
    df : polars.DataFrame
        The input data containing processed variant information with chromosome and density details.

    sorted_value : list
        A list of chromosome labels to specify the order in which the x-axis should be sorted.

    save_name : str, optional (default: 'Figure1B')
        The filename used to save the generated plot. The plot will be saved in both .svg and .png formats.

    y : str, optional (default: 'Percent_den:Q')
        The column to be used for the y-axis, representing the percentage density of variants.

    x : str, optional (default: 'column_3_right')
        The column to be used for the x-axis, representing the chromosome labels.

    color : str, optional (default: 'column_3_right:N')
        The column to be used for color encoding, representing the chromosomes.

    title : str, optional (default: 'Variant occurrence per chromosome with resolution of 10Kb')
        The title of the plot.

    Returns:
    --------
    altair.Chart
        The generated Altair scatter plot.

    Notes:
    ------
    - The plot displays variant occurrences with a jitter effect to avoid overlapping points.
    - The generated plot uses category20b colors and is saved at a resolution suitable for publication.

    Example Usage:
    --------------
    >>> plot = Plot_VarChr(df, sorted_value=['1', '2', '3', 'X', 'Y'], save_name='ChromosomeVariantPlot')
    
    """
    
    # Set color and sort labels
    color_scale = alt.Scale(scheme="category20b")

    # Create plot
    plot = alt.Chart(df, title=title).mark_circle(size=8).encode(
        y=alt.Y(y, title="Percentage of Variant"),
        x=alt.X(x, title="Chromosome", axis=alt.Axis(grid=True, labelAngle=-45, ticks=False),
                sort=sorted_value),
        xOffset="jitter:Q",
        color=alt.Color(color, legend=None, scale=color_scale)
    ).transform_calculate(
        jitter="sqrt(-2*log(random()))*cos(2*PI*random())"
    )

    # Save
    plot.save(f'{save_name}.svg', engine='vl-convert', ppi=300)
    plot.save(f'{save_name}.png', engine='vl-convert', ppi=300)
    
    return plot


