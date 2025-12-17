"""
Constructs the working cohort by doing the following:
- Load the ClinVar file (variant_summary.txt)
- Filter for GRCh38 (current human genome)
- Filter for Missense variants
- Filter for Pathogenic or Benign
- Map the ClinVar Gene Symbols to UniProt IDs
- Save a clean cohort.csv
"""

import pandas as pd
import os
import re
from datetime import datetime

# CONFIGURATION #
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLINVAR_PATH = os.path.join(BASE_DIR, 'data', 'raw', 'variant_summary.txt')
MAPPING_PATH = os.path.join(BASE_DIR, 'data', 'raw', 'human_id_mapping.tsv')
OUTPUT_PATH = os.path.join(BASE_DIR, 'data', 'processed', 'cohort_mapped.csv')

# Create processed folder if it doesn't exist
os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

def parse_amino_acid_change(name):
    """
    Extracts the amino acid (AA) change from the ClinVar Name string.
    Expected format: ... (p.Glu123Asp) ...
    Returns: (WildType, ResidueIndex, Mutant) or (None, None, None)
    """
    # Regex search method:
    # p\.               -> Looks for preceding p.
    # ([A-Z][a-z]{2})   -> Matches 3-letter AA code
    # (\d+)             -> Matches residue number
    # ([A-Z][a-z]{2})   -> Matches 3-letter AA code
    match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', str(name))
    if match:
        return match.group(1), int(match.group(2)), match.group(3)
    return None, None, None

def build_cohort():
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Loading ClinVar data...")

    df = pd.read_csv(CLINVAR_PATH, sep='\t', low_memory=False)
    print(f"{datetime.now().strftime('%H:%M:%S')}: ClinVar data loaded.")

    initial_count = len(df)
    print(f"\nInitial ClinVar rows: {initial_count}")

    # FILTER 1: GENOME ASSEMBLY #
    # We only want the latest human genome build
    df = df[df['Assembly'] == 'GRCh38']
    print(f"Rows after Assembly 'GRCh38' filter: {len(df)}")
    
    # FILTER 2: VARIANT TYPE #
    # We only want Single Nucleotide Variants (SNVs)
    df = df[df['Type'] == 'single nucleotide variant']
    print(f"Rows after SNV filter: {len(df)}")

    # FILTER 3: REVIEW STATUS #
    # We reject 'no assertion criteria provided' or 'no assertion for the individual variant'
    # We accept: 'criteria provided...', 'reviewed by expert panel', 'practice guideline'
    bad_status = ['no assertion criteria provided', 'no assertion for the individual variant']
    df = df[~df['ReviewStatus'].isin(bad_status)]
    print(f"Rows after ReviewStatus filter: {len(df)}")

    # FILTER 4: MISSENSE EXTRACTION #
    # Apply the regex function to the 'Name' column, creating a DataFrame with 3 new columns
    aa_data = df['Name'].apply(lambda x: pd.Series(parse_amino_acid_change(x)))
    aa_data.columns = ['WildType', 'ResidueIndex', 'MutantAA']

    # Join these new columns to our main DataFrame
    df = pd.concat([df, aa_data], axis=1)

    # Drop rows where Regex failed (meaning it wasn't a standard missense mutation)
    df = df.dropna(subset=['WildType', 'ResidueIndex'])
    print(f"Rows after Missense Extraction: {len(df)}")

    # FILTER 5: CLINICAL SIGNIFICANCE #
    # Group into 0 (Benign) and 1 (Pathogenic)
    # Ignore VUS (Uncertain significance)

    def simplify_significance(sig):
        sig = str(sig).lower()
        if 'pathogenic' in sig and 'benign' not in 'sig' and 'conflict' not in sig:
            return 'Pathogenic'
        if 'benign' in sig and 'pathogenic' not in 'sig' and 'conflict' not in sig:
            return 'Benign'
        return 'Other'
    
    df['Class'] = df['ClinicalSignificance'].apply(simplify_significance)
    df = df[df['Class'] != 'Other']
    print(f"Rows after Class filter (Pathogenic/Benign): {len(df)}")

    # FILTER 6: MAPPING TO UNIPROT #
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Loading UniProt Mapping file...")

    # Loading just the needed columns: 'Gene Names' and 'Entry' (UniProt ID)
    # UniProt columns are: Entry, Reviewed, Entry Name, Protein names, Gene Names, Organism, Length, PDB, AlphaFoldDB, RefSeq
    try:
        mapping_df = pd.read_csv(MAPPING_PATH, sep='\t')

        # Need to map ClinVar's 'GeneSymbol' to UniProt 'Gene Names'
        # UniProt 'Gene Names' look like "TP53 P53", but we will do a merge

        # Explode UniProt genes (if multiple names exist) or just take the first word
        mapping_df['GeneSymbol'] = mapping_df['Gene Names'].str.split(' ').str[0]

        # Keep only the ID and Gene, then rename
        mapping_df = mapping_df[['Entry', 'GeneSymbol']]
        mapping_df = mapping_df.rename(columns={'Entry': 'UniProtID'})

        # Merge
        merged_df = pd.merge(df, mapping_df, on='GeneSymbol', how='inner')
        print(f"Rows after Mapping to UniProt IDs: {len(merged_df)}")

    except Exception as e:
        print(f"WARNING: Mapping failed: {e}")
        print("Saving unmapped cohort instead.")
        merged_df = df

    # SAVE #
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Saving final cohort to {OUTPUT_PATH}")
    columns_to_keep = ['Name', 'GeneSymbol', 'UniProtID', 'Chromosome', 'WildType', 'ResidueIndex', 'MutantAA', 'Class', 'ReviewStatus']
    final_cols = [c for c in columns_to_keep if c in merged_df.columns]

    merged_df[final_cols].to_csv(OUTPUT_PATH, index=False)
    print("Saved final cohort.")

if __name__ == "__main__":
    build_cohort()