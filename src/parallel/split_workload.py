"""
Splits workload by doing the following:
- Loads cohort_filtered.csv
(- We know that repairing the protein is the most intensive step, so we should repair the proteins with the most mutations first)
- Sorts cohort by efficiency (# of mutations for protein)
- Cuts into chunks for parallel processing
(- This computer has 8 cores, 4 for performance, 4 for efficiency. For safety, we'll use 6/8)
"""

import pandas as pd
import os
import numpy as np

# CONFIGURATION #
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# INPUT_CSV = os.path.join(BASE_DIR, 'data', 'processed', 'cohort_filtered.csv')
INPUT_CSV = os.path.join(BASE_DIR, 'data', 'processed', 'cohort_part_6.csv')
NUM_CHUNKS = 6

def split_csv():
    print(f"Loading {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)

    # Count mutations per protein
    counts = df['UniProtID'].value_counts()
    # Map counts back to the dataframe
    df['MutationCount'] = df['UniProtID'].map(counts)

    # Sort descending: proteins with more mutations go first
    df_sorted = df.sort_values(by=['MutationCount', 'UniProtID'], ascending=[False, True])
    df_sorted = df_sorted.drop(columns=['MutationCount']) # Drop the helper column
    
    print(f"Sorted {len(df_sorted)} rows. High-yield proteins are now first.")

    # Chop dataframe into N parts
    chunks = np.array_split(df_sorted, NUM_CHUNKS)
    
    for i, chunk in enumerate(chunks):
        worker_id = i + 1
        output_name = f"cohort_part_6_{worker_id}.csv"
        output_path = os.path.join(BASE_DIR, 'data', 'processed', output_name)
        
        chunk.to_csv(output_path, index=False)
        print(f"  > Created {output_name} ({len(chunk)} rows) - For Worker {worker_id}")

if __name__ == "__main__":
    split_csv()