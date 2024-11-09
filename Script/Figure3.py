import altair as alt
alt.data_transformers.enable("vegafusion")
import polars as pl

def plot_sv_chart(data_path, save_name='Figure3BC'):
    """
    Generates combined bar charts for Median Length and Total Length of Structural Variants.

    The function reads a TSV file containing data about structural variants and creates two bar charts: 
    one for the median length and one for the total length of the structural variants. These charts are 
    then combined and saved as SVG and PNG files.

    Parameters:
    data_path (str): Path to the input TSV data file. The file should contain columns 'SVTYPE', 'Metric', and 'Value'.
    save_name (str): The base name to save the output SVG and PNG charts.

    Returns: combined_chart
    altair.Chart
        The generated Altair bar chart.
    """
    # Read data
    med_len_SV = pl.read_csv(data_path, separator='\t')

    # Base chart configuration
    base = alt.Chart(med_len_SV, title=" ").encode(
        alt.X('SVTYPE:N', title='SV Type', axis=alt.Axis(labelAngle=0, grid=True)),
        color=alt.Color('Metric:N', title=None, scale=alt.Scale(scheme="set1"), sort=["Median Length (bp)","Total Length (Mbp)"])
    )

    # Median Length bar chart
    median_bar = base.transform_filter(
        alt.datum.Metric == 'Median Length (bp)'
    ).mark_bar(opacity=0.7).encode(
        alt.Y('Value:Q', axis=alt.Axis(title='Median Length (bp)')),
        xOffset='Metric:N'
    )

    # Total Length bar chart
    total_bar = base.transform_filter(
        alt.datum.Metric == 'Total Length (Mbp)'
    ).mark_bar(opacity=0.7).encode(
        alt.Y('Value:Q', axis=alt.Axis(title='Total Length (Mbp)', offset=0)),
        xOffset='Metric:N'
    )

    # Combine the bar charts
    combined_chart = median_bar | total_bar

    # Save charts
    combined_chart.save(f"{save_name}.svg", engine="vl-convert", ppi=300)
    combined_chart.save(f"{save_name}.png", engine="vl-convert", ppi=300)
    
    return combined_chart


