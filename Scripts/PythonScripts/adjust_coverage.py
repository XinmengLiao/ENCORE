import pandas as pd

# Read the coverage table
df = pd.read_csv('coverage_table.tsv', sep='\t')

# Remove .concoct_part_0 from first column (skip header)
df.iloc[1:, 0] = df.iloc[1:, 0].str.replace('.concoct_part_0', '')

# Save modified table
df.to_csv('coverage_table.tsv', sep='\t', index=False)