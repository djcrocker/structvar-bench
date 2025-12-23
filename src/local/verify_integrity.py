"""
Checks proper file and column loading for ClinVar and AlphaFold (1 sample) files.
"""

import pandas as pd
import gzip
from Bio.PDB import MMCIFParser, PDBParser
import os

# CONFIGURATION #
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLINVAR_PATH = os.path.join(BASE_DIR, 'data', 'raw', 'variant_summary.txt')
ALPHAFOLD_TEST_FILE = os.path.join(BASE_DIR, 'data', 'raw', 'alphaFold_human', 'AF-A0A0A0MRZ7-F1-model_v6.cif.gz')

def check_clinvar():
    print(f"Checking ClinVar Data: {CLINVAR_PATH}")

    if not os.path.exists(CLINVAR_PATH):
        print(f"ERROR: File not found at {CLINVAR_PATH}")
        return
    
    try:
        df = pd.read_csv(CLINVAR_PATH, sep='\t', low_memory=False) # Tab-separated, low_memory=False handles DtypeWarning
        print(f"SUCCESS: File loaded successfully. Total Rows: {len(df)}")
        
        # Check for required columns
        required_cols = ['Name', 'Chromosome', 'ClinicalSignificance', 'ReviewStatus']
        missing = [col for col in required_cols if col not in df.columns]
        
        if missing:
            print(f"FAILURE: Missing columns: {missing}")
        else:
            print(f"SUCCESS: All required columns found: {required_cols}")

        print("\nSnapshot of first 5 variants:")
        print(df[required_cols].head())
        
    except Exception as e:
        print(f"FAILURE: Failed to parse ClinVar file: {e}")

def check_alphafold():
    print(f"Checking AlphaFold Structure with test file: {ALPHAFOLD_TEST_FILE}")

    if not os.path.exists(ALPHAFOLD_TEST_FILE):
        print(f"ERROR: Test file not found at {ALPHAFOLD_TEST_FILE}")
        return
    
    try:
        if ALPHAFOLD_TEST_FILE.endswith('.cif.gz'):
            parser = MMCIFParser(QUIET=True)
        else: # .pdb.gz
            parser = PDBParser(QUIET=True)

        # Parse the structure using MMCIFParser/PDBParser + gzip
        # Open the .gz file in 'rt' (Read Text) mode and pass the handle to Biopython
        with gzip.open(ALPHAFOLD_TEST_FILE, 'rt') as handle:
            # Bioython reads directly from the unzipped stream
            structure = parser.get_structure('test_protein', handle)
        print("SUCCESS: .gz file parsed successfully.")
        
        # Check for pLDDT Scores (B-factors)
        first_model = structure[0]
        first_chain = list(first_model)[0]
        first_residue = list(first_chain)[0]
        first_atom = list(first_residue)[0]
        
        b_factor = first_atom.bfactor
        
        print(f"SUCCESS: B-factor (pLDDT) found. Value for first atom: {b_factor}")
        
        if 0 <= b_factor <= 100:
            print("VALIDATION: pLDDT score is within expected range (0-100).")
        else:
            print(f"WARNING: Unusual pLDDT score: {b_factor}. Check if this is a standard AF output.")

    except Exception as e:
        print(f"FAILURE: Failed to parse file: {e}")

if __name__ == "__main__":
    print("Starting Integrity Verification\n")
    check_clinvar()
    check_alphafold()
    print("\nVerification Complete")