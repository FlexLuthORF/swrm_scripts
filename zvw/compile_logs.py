import pandas as pd
import re

# Function to extract CONSCOUNT and isotype from FIELD2 in assemble_table
def extract_assemble_data(compound_field):
    match = re.search(r"CONSCOUNT=(\d+)\|PRCONS=([A-Z0-9]+)", compound_field)
    if match:
        return int(match.group(1)), match.group(2)
    return 0, None

# Function to extract CONSCOUNT and isotype from ID in maskquals_table
def extract_maskquals_data(compound_field):
    match = re.search(r"\|CONSCOUNT=([0-9,]+)\|PRCONS=([A-Z0-9]+)", compound_field)
    if match:
        conscount = sum(map(int, match.group(1).split(',')))
        return conscount, match.group(2)
    return 0, None

# Initialize an empty DataFrame to store results
result_df = pd.DataFrame()

# List of files to process
files = ["assemble_table.tab", "cregion_table.tab", "consensus-crr_table.tab", 
         "consensus-vrr_table.tab", "maskqual_table.tab", "primers-vrr_table.tab", 
         "final-unique_headers.tab", "final-total_headers.tab", "final-unique-atleast2_headers.tab"]

# Process each file
for file_name in files:
    df = pd.read_csv(file_name, sep='\t', low_memory=False)

    if file_name in ["cregion_table.tab", "primers-vrr_table.tab"]:
        # Count rows per SampleID
        counts = df.groupby('SampleID').size() - 1
        result_df = pd.concat([result_df, counts.rename(file_name)], axis=1)

    elif file_name in ["consensus-vrr_table.tab", "consensus-crr_table.tab"]:
        # Sum CONSCOUNT per SampleID and isotype
        aggregated = df.groupby(['SampleID', 'PRCONS'])['CONSCOUNT'].sum().unstack(fill_value=0)
        result_df = pd.concat([result_df, aggregated.add_suffix('_' + file_name)], axis=1)

    elif file_name == "assemble_table.tab":
        if 'FIELD2' in df.columns:
            df['CONSCOUNT'], df['Isotype'] = zip(*df['FIELD2'].apply(extract_assemble_data))
            aggregated = df.groupby(['SampleID', 'Isotype'])['CONSCOUNT'].sum().unstack(fill_value=0)
            result_df = pd.concat([result_df, aggregated.add_suffix('_' + file_name)], axis=1)

    elif file_name == "maskqual_table.tab":
        if 'ID' in df.columns:
            df['CONSCOUNT'], df['Isotype'] = zip(*df['ID'].apply(extract_maskquals_data))
            aggregated = df.groupby(['SampleID', 'Isotype'])['CONSCOUNT'].sum().unstack(fill_value=0)
            result_df = pd.concat([result_df, aggregated.add_suffix('_' + file_name)], axis=1)

    elif file_name.startswith("final"):
        # Use CONSCOUNT and PRCONS directly
        aggregated = df.groupby(['SampleID', 'PRCONS'])['CONSCOUNT'].sum().unstack(fill_value=0)
        result_df = pd.concat([result_df, aggregated.add_suffix('_' + file_name)], axis=1)

# Fill NaNs with zeros and reset index
result_df = result_df.fillna(0).reset_index()

# Save to a new file
result_df.to_csv("compiled_results.csv", index=False)
