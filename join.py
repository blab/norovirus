# Import pandas library
import pandas as pd

# Read first two csv files with '\t' separator
tsv1 = pd.read_csv("config/reference_sequences.tsv", sep='\t')
tsv2 = pd.read_csv("results/metadata_parsed_big.tsv", sep='\t')

# store the result in Output_df dataframe.
# Here common column is 'ID' column
Output_df = pd.merge(tsv1, tsv2, on='strain',
                     how='inner')

# Now store the 'Output_df'
# in tsv file 'Output.tsv'
Output_df.to_csv("results/metadata_parsed.tsv",
                 sep="\t", header=True,
                 index=False)
