"""
Cross-references structures from the mapped cohort with the AlphaFold DB.
"""

import pandas as pd
import os
import gzip
from Bio.PDB import MMCIFParser, PDBParser
import warnings
from datetime import datetime

# CONFIGURATION #
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INPUT_CSV = os.path.join(BASE_DIR, 'data', 'processed', 'cohort_mapped.csv')
OUTPUT_CSV = os.path.join(BASE_DIR, 'data', 'processed', 'cohort_filtered.csv')
ALPHAFOLD_DIR = os.path.join(BASE_DIR, 'data', 'raw', 'alphaFold_human')
MIN_PLDDT = 70.0

def get_structure_plddt(filepath, residue_index):
    """
    Picks the right file parser based on file type and returns the pLDDT.
    """
    try:
        # Choose the Parser
        if filepath.endswith('.cif') or filepath.endswith('.cif.gz'):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        # Handle gzip if needed
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rt') as handle:
                structure = parser.get_structure('temp', handle)
        else:
            structure = parser.get_structure('temp', filepath)

        # Get the model/chain
        model = structure[0]
        chain = list(model)[0]

        # Find Residue
        if residue_index in chain:
            residue = chain[residue_index]
            if 'CA' in residue:
                return residue['CA'].bfactor
            else:
                # Fallback = average of all atoms
                bfactors = [atom.bfactor for atom in residue]
                return sum(bfactors) / len(bfactors)
        else:
            return None 
            
    except Exception as e:
        print(f"ERROR: Error parsing {filepath}: {e}")
        return None
    
def filter_cohort():
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Loading {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)

    print(f"{datetime.now().strftime('%H:%M:%S')}: Scanning {ALPHAFOLD_DIR} for structure files...")
    available_files = os.listdir(ALPHAFOLD_DIR)
    print(f"{datetime.now().strftime('%H:%M:%S')}: Found {len(available_files)} available structure files.")

    # Map UniProtID -> Filename (Handling PDB and CIF)
    file_map = {}
    for f in available_files:
        # Looks for files starting with AF- and having a valid extension
        if f.startswith('AF-') and ('model' in f): 
            # Parse ID (AF-P12345-F1...)
            parts = f.split('-')
            if len(parts) >= 2:
                uniprot_id = parts[1]
                # If we have a duplicate (e.g. both .cif and .pdb), prefer .pdb for FoldX compatibility
                if uniprot_id not in file_map:
                    file_map[uniprot_id] = f
                elif f.endswith('.pdb') or f.endswith('.pdb.gz'): # if ID already in map, only overwrite if new one is .pdb
                    file_map[uniprot_id] = f

    print(f"Found {len(file_map)} unique matching structure files.")

    df['StructureFile'] = df['UniProtID'].map(file_map)

    # MISSING FILES - Step 2.II #
    before_count = len(df)
    df = df.dropna(subset=['StructureFile'])
    print(f"Dropped {before_count - len(df)} rows missing structure files.")

    # LOW pLDDT SCORE - Step 2.III#

    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Checking pLDDT scores...")
    
    plddt_scores = []
    
    total = len(df)
    for i, (index, row) in enumerate(df.iterrows()):
        if i % 100 == 0: print(f"Processing {i}/{total}...", end='\r')
        
        full_path = os.path.join(ALPHAFOLD_DIR, row['StructureFile'])
        score = get_structure_plddt(full_path, row['ResidueIndex'])
        plddt_scores.append(score)

    df['pLDDT'] = plddt_scores
    
    # Drop Low Confidence
    df_clean = df[df['pLDDT'] >= MIN_PLDDT].copy()
    
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Processing Complete.")
    print(f"Final Cohort Size: {len(df_clean)}")

    df_clean.to_csv(OUTPUT_CSV, index=False)
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Saving final cohort to {OUTPUT_CSV}")

if __name__ == "__main__":
    filter_cohort()