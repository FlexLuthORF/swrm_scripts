import matplotlib.pyplot as plt

# Pivoting the data to get the sum of CONSCOUNT for each PRCONS per SampleID
pivoted_data = data.pivot_table(index='SampleID', columns='PRCONS', values='CONSCOUNT', aggfunc='sum', fill_value=0)

# Color mapping provided by the user
color_mapping = {
    'IgG': '#0000FF',  # Blue
    'IgM': '#008000',  # Green
    'IgK': '#FFA500',  # Orange
    'IgL': '#800080',  # Purple
    'IgA': '#FF0000',  # Red
    'IgD': '#00CED1',  # DarkTurquoise
}

# Define the number of pie charts per row
charts_per_row = 4

# Function to generate pie charts
def generate_pie_charts(pivoted_data, color_mapping, charts_per_row):
    # Number of rows and columns for the subplot grid
    num_samples = len(pivoted_data)
    num_rows = -(-num_samples // charts_per_row)  # Ceiling division
    num_cols = charts_per_row

    # Creating figure with subplots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(4 * num_cols, 4 * num_rows))
    axs = axs.flatten()  # Flattening the axis array

    # Hide axes for subplots without a pie chart
    for i in range(num_samples, len(axs)):
        axs[i].axis('off')

    # Generating each pie chart
    for i, (sample_id, row) in enumerate(pivoted_data.iterrows()):
        # Colors for each PRCONS
        colors = [color_mapping.get(prcons, 'grey') for prcons in row.index]

        # Calculating percentages and labels for pie slices
        total_count = row.sum()
        percentages = row / total_count * 100
        labels = [f'{prcons} ({count})' for prcons, count in zip(row.index, row)]

        # Plotting pie chart
        axs[i].pie(percentages, labels=None, colors=colors, autopct='%1.1f%%', startangle=140)
        axs[i].set_title(sample_id)

    # Adjusting layout
    plt.tight_layout()

    return fig, axs

# Generating the pie charts
fig, axs = generate_pie_charts(pivoted_data, color_mapping, charts_per_row)

# Adding a legend outside of the last ax
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right')

# Saving the figure
output_file = '/mnt/data/pie_charts_custom.png'
fig.savefig(output_file, bbox_inches='tight', dpi=300)

output_file
