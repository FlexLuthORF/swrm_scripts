import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = 'cov_merged.tsv'
df = pd.read_csv(file_path, sep='\t', header=None)

# Name columns for easier reference
df.columns = ['region', 'column2', 'coverage', 'column4', 'column5', 'column6', 'sampleId']

# Function to determine the grouping of regions
def group_region(region):
    if region.startswith('IGH'):
        if len(region) > 3 and region[3] in 'VDJ':
            return 'IGH'
        else:
            return 'IGHC'
    return region[:3]

# Apply grouping function to regions
df['grouped_region'] = df['region'].apply(group_region)

# Compute the average coverage for each sampleId and grouped region
avg_coverage = df.groupby(['sampleId', 'grouped_region'])['coverage'].mean().reset_index()

# Save the computed averages to a .tsv file
output_path = 'avg_coverage.tsv'
avg_coverage.to_csv(output_path, sep='\t', index=False)

# Pivot the data for plotting
pivot_df = avg_coverage.pivot(index='grouped_region', columns='sampleId', values='coverage')

# Plot the data
pivot_df.plot(kind='bar', figsize=(12, 8))
plt.title('Average Coverage by Region and Sample')
plt.xlabel('Region')
plt.ylabel('Average Coverage')
plt.axhline(y=50, color='r', linestyle='--')
plt.axhline(y=20, color='b', linestyle='--')
plt.legend().set_visible(False)
plt.tight_layout()

# Save plot
plot_path_with_lines = 'avg_coverage_plot_with_lines.png'
plt.savefig(plot_path_with_lines)
