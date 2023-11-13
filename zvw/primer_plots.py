import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the table
cregion_table = pd.read_csv("logs/cregion_table.tab", delimiter="\t")

# Print frequency table of the "PRIMER" column
print(cregion_table['c_call'].value_counts())

# Create a color palette
color_palette = {
    "IGHA": "#882255",
    "IGHD": "#AA4499",
    "IGHE": "#88CCEE",
    "IGHG": "#CC6677",
    "IGHM": "#6699CC",
    "IGKC": "#44AA99",
    "IGLC1": "#888888",
    "IGLC3": "#DDCC77"
}

# Sort and assign the names to the color palette
isotypes = sorted(cregion_table['PRIMER'].unique())
color_palette = {key: color_palette[key] for key in isotypes}

# Plot isotype primer position
sns.set_theme(style="whitegrid")
plt.figure(figsize=(12, 6))
sns.histplot(
    cregion_table,
    x="PRSTART",
    hue="PRIMER",
    element="step",
    bins=1,
    palette=color_palette
)
plt.xlabel("C-region alignment start")
plt.ylabel("Count")
plt.title("Isotype Primer Position")
plt.legend(title="Primer", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('primer_plot.png')
plt.show()

