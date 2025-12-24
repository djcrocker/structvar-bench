"""
Performs structure repair and energy extraction by doing the following (call the protein A12345):
- Converts A12345.cif.gz or A12345.pdb.gz to A12345.pdb for FoldX
- Runs FoldX RepairPDB (thermodynamic minimization) on A12345.pdb, resulting in A12345_Repair.pdb (saves repaired A12345-OriginalResidueMutant.pdb to /data/processed/structures/)
- Runs FoldX BuildModel (mutation) on A12345_Repair.pdb, resulting in Dif_A12345_Repair.fxout
- Looks through .fxout for ddG value, adds value to cohort_with_ddg.csv
Usage: python3 src/parallel/worker_foldx.py --worker_id # --batch_size #
"""

import pandas as pd
import os
import subprocess
import gzip
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
import shutil
from datetime import datetime
import argparse

# ARGUMENT PARSING #
# Lets us run: python3 src/parallel/worker_foldx.py --worker_id 1
parser = argparse.ArgumentParser()
parser.add_argument('--worker_id', type=str, default='1', help='ID for this parallel worker (1, 2, 3...)')
parser.add_argument('--batch_size', type=int, default='1', help='Number of rows to process')
args = parser.parse_args()
WORKER_ID = args.worker_id
BATCH_LIMIT = args.batch_size

# GLOBAL CONFIGURATION
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ALPHAFOLD_DIR = os.path.join(BASE_DIR, 'data', 'raw', 'alphaFold_human') 

# INPUT/OUTPUT (depends on the Worker ID) #
# e.g. cohort_part_1.csv -> cohort_with_ddg_1.csv
INPUT_CSV = os.path.join(BASE_DIR, 'data', 'processed', f'cohort_part_{WORKER_ID}.csv')
OUTPUT_CSV = os.path.join(BASE_DIR, 'data', 'processed', f'cohort_with_ddg_{WORKER_ID}.csv')

# FOLDX PATH #
FOLDX_BIN = os.path.join(BASE_DIR, 'tools', 'foldx', 'foldx')

# UNIQUE WORKER WORKSPACES #
WORK_DIR = os.path.join(BASE_DIR, 'data', 'processed', f'foldx_workspace_{WORKER_ID}')
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, 'data', 'processed', 'structures'), exist_ok=True)

# PILOT MODE #
PILOT_MODE = False
PILOT_SIZE = 5

class ChainASelect(Select):
    """AlphaFold models often have confident regions in Chain A. 
    Filter to ensure only Chain A is saved to keep FoldX happy."""
    def accept_chain(self, chain):
        return chain.id == 'A'
    
def convert_cif_to_pdb(structure_path, output_pdb_path):
    """Converts .cif.gz or .pdb.gz to .pdb for FoldX."""
    try:
        # Select the right parser
        if '.pdb' in structure_path:
            parser = PDBParser(QUIET=True)
        else:
            parser = MMCIFParser(QUIET=True)

        # Opening file by handling GZIP
        if structure_path.endswith('.gz'):
            with gzip.open(structure_path, 'rt') as f:
                structure = parser.get_structure('temp', f)
        else:
            structure = parser.get_structure('temp', structure_path)
        
        # Now save as clean PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_path, select=ChainASelect())
        return True
    except Exception as e:
        print(f"Error converting {structure_path}: {e}")
        return False
    
def get_1letter_code(aa3: str):
    """Converts 3-letter AA code to 1-letter code. Handles exceptions."""
    return seq1(aa3)
    
def run_foldx_process(df):
    # Load existing results
    completed_keys = set()
    write_header = True
    if os.path.exists(OUTPUT_CSV):
        # Read existing file to see what has already been finished
        try:
            print(f"\n{datetime.now().strftime('%H:%M:%S')}: Found existing output file. Scanning for completed variants.")

            # Ensure that new results aren't appended to middle of last line
            with open(OUTPUT_CSV, 'rb+') as f:
                f.seek(0, 2)  # Move to the very end of the file
                size = f.tell()
                if size > 0:
                    f.seek(-1, 2)  # Move back 1 character
                    last_char = f.read(1)
                    if last_char != b'\n':
                        print("  > Fixed missing newline at end of CSV.")
                        f.write(b'\n') # Force add the newline

            existing_df = pd.read_csv(OUTPUT_CSV)
            # Create a unique key for every finished row: UniprotID_WildTypeResNumMutantAA
            for _, r in existing_df.iterrows():
                key = f"{r['UniProtID']}_{r['WildType']}{int(r['ResidueIndex'])}{r['MutantAA']}"
                completed_keys.add(key)
            print(f"Resuming: {len(completed_keys)} variants already finished.")
            write_header = False # File exists, so don't write header again
        except Exception as e:
            print(f"ERROR: Couldn't read CSV ({e}). Starting fresh.")
    
    # Group by Protein (UniProtID) to minimize Repair steps
    grouped = df.groupby('UniProtID')
    total_proteins = len(grouped)

    variants_processed_session = 0
    
    print(f"\n{datetime.now().strftime('%H:%M:%S')}: Starting processing...")
    if BATCH_LIMIT:
        print(f"Batch Limit Active: Will stop after {BATCH_LIMIT} new mutations.")
    
    for i, (uniprot_id, group) in enumerate(grouped):
        # Check if protein is fully done
        group_keys = [f"{uniprot_id}_{row['WildType']}{int(row['ResidueIndex'])}{row['MutantAA']}" for _, row in group.iterrows()]
        if all(k in completed_keys for k in group_keys):
            # print(f"  > Skipping {uniprot_id} (All variants complete).")
            continue
        
        print(f"\n{datetime.now().strftime('%H:%M:%S')}: [{i+1}/{total_proteins}] Processing {uniprot_id} ({len(group)} mutations)...")

        # PREPARE THE PDB - Step 3.II #
        # Look for the source structure
        struct_file = group.iloc[0]['StructureFile']
        cif_path = os.path.join(ALPHAFOLD_DIR, struct_file)
        
        # Define local paths in workspace
        raw_pdb_name = f"{uniprot_id}.pdb"
        repaired_pdb_name = f"{uniprot_id}_Repair.pdb" # FoldX naming convention
        
        # Check if already repaired
        if not os.path.exists(os.path.join(WORK_DIR, repaired_pdb_name)):
            # Convert CIF -> PDB
            print(f"  > Converting CIF to PDB...")
            if not convert_cif_to_pdb(cif_path, os.path.join(WORK_DIR, raw_pdb_name)):
                print("  > Conversion failed. Skipping protein.")
                continue

            # Run FoldX Repair (thermodynamic minimization) <- Step 3.II
            print(f"  > Running FoldX Repair...")
            cmd = [
                FOLDX_BIN, 
                '--command=RepairPDB', 
                f'--pdb={raw_pdb_name}',
                '--ionStrength=0.05',
                '--pH=7.0',
                '--vdwDesign=2'
            ]
            
            # Run inside WORK_DIR so output files stay there
            try:
                subprocess.run(cmd, cwd=WORK_DIR, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                # subprocess.run(cmd, cwd=WORK_DIR, check=True) # Use to see FoldX working
            except subprocess.CalledProcessError:
                print("  > FoldX Repair crashed. Skipping.")
                continue
        else:
            print("  > Repaired structure found. Using cached version.")

        print(f"\n{datetime.now().strftime('%H:%M:%S')}: Protein repaired. Running mutations...")

        # RUN MUTATIONS - Step 3.III #
        MUT_LIST_FILE = "individual_list.txt"
        results_batch = []
        for idx, row in group.iterrows():
            # Build unique key again to check for variants
            res_num = int(row['ResidueIndex'])
            unique_key = f"{uniprot_id}_{row['WildType']}{res_num}{row['MutantAA']}"

            if unique_key in completed_keys:
                continue

            wt_1 = get_1letter_code(row['WildType'])
            mut_1 = get_1letter_code(row['MutantAA'])
            
            # Construct Mutation String (e.g. "EA10D" - OriginalAA Chain ResNum TargetAA)
            # AlphaFold is almost always Chain A
            mut_code = f"{wt_1}A{res_num}{mut_1}"
            
            # Config file for FoldX
            with open(os.path.join(WORK_DIR, MUT_LIST_FILE), 'w') as f:
                f.write(f"{mut_code};")

            print(f"  > Mutating {mut_code}...", end="")
            
            cmd = [
                FOLDX_BIN,
                '--command=BuildModel',
                f'--pdb={repaired_pdb_name}',
                f'--mutant-file={MUT_LIST_FILE}',
                '--numberOfRuns=1' # For publications, usually 3-5. Here, keep it fast.
            ]
            
            try:
                subprocess.run(cmd, cwd=WORK_DIR, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                # subprocess.run(cmd, cwd=WORK_DIR, check=True) # Use to see FoldX working
                
                # PARSE ENERGIES - Step 3.IV #
                # FoldX appends results to Dif_[PDBname]_Repair.fxout
                fxout_name = f"Dif_{repaired_pdb_name.replace('.pdb', '.fxout')}"
                fxout_path = os.path.join(WORK_DIR, fxout_name)

                ddg_val = None
                
                if os.path.exists(fxout_path):
                    with open(fxout_path, 'r') as f:
                        lines = f.readlines()
                        if len(lines) > 0:
                            last_line = lines[-1]
                            cols = last_line.split('\t')
                            # Ensure it's a data line, not a header
                            if len(cols) > 2 and "Total" not in cols[1]:
                                try:
                                    ddg_val = float(cols[1])
                                except ValueError:
                                    pass

                # Rename and organize mutant structure
                # FoldX generates "P12345_Repair_1.pdb", we want P12345_A10D.pdb
                generated_pdb = os.path.join(WORK_DIR, f"{repaired_pdb_name.replace('.pdb', '_1.pdb')}")
                final_pdb_name = f"{uniprot_id}_{wt_1}{res_num}{mut_1}.pdb"
                final_pdb_path = os.path.join(BASE_DIR, 'data', 'processed', 'structures', final_pdb_name)

                os.makedirs(os.path.dirname(final_pdb_path), exist_ok=True)
                
                if os.path.exists(generated_pdb):
                    shutil.move(generated_pdb, final_pdb_path)
                else:
                    print(f"  > WARNING: Mutant PDB not generated for {mut_code}")

                # Record results
                result_entry = row.to_dict()
                result_entry['ddG'] = ddg_val
                result_entry['MutantStructureFile'] = final_pdb_name
                results_batch.append(result_entry)

                completed_keys.add(unique_key)
                variants_processed_session += 1

                print(f"\n  > Done. ddG: {ddg_val}")
                print(f"Processed variants in session: {variants_processed_session}")

                # Clean up intermediate files
                try:
                    # Files generated by BuildModel
                    trash = [
                        fxout_name,                                 # Dif_...fxout
                        fxout_name.replace('Dif_', 'Raw_'),         # Raw_...fxout
                        fxout_name.replace('Dif_', 'PdbList'),      # PdbList_...fxout
                        fxout_name.replace('Dif_', 'Average'),      # Average_...fxout
                        'individual_list.txt'
                    ]

                    # Also look for "WT_" artifact that FoldX sometimes makes
                    wt_artifact = os.path.join(WORK_DIR, f"WT_{repaired_pdb_name.replace('.pdb', '_1.pdb')}")
                    if os.path.exists(wt_artifact):
                        trash.append(os.path.basename(wt_artifact))

                    for t in trash:
                        full_path = os.path.join(WORK_DIR, t)
                        if os.path.exists(full_path):
                            os.remove(full_path)
                            
                except Exception as e:
                    print(f"Cleanup warning: {e}")

            except subprocess.CalledProcessError:
                print(f"\n  > FoldX BuildModel failed for {mut_code}. Skipping.")

        # Save results after every protein
        if results_batch:
            batch_df = pd.DataFrame(results_batch)
            batch_df.to_csv(OUTPUT_CSV, mode='a', header=write_header, index=False)
            
            # CRITICAL: After writing the first batch, turn off header for future batches
            write_header = False 
            print(f"BATCH SAVED ({len(results_batch)} mutations). Session Total: {variants_processed_session}")

        if BATCH_LIMIT and variants_processed_session >= BATCH_LIMIT:
            print(f"\n{datetime.now().strftime('%H:%M:%S')}: Batch limit of {BATCH_LIMIT} reached (Processed {variants_processed_session}). Stopping safely.")
            break

if __name__ == "__main__":
    # Load Data
    df = pd.read_csv(INPUT_CSV)
    print("--- FoldX Pipeline ---")

    if PILOT_MODE:
        print(f"Pilot Mode Active: First {PILOT_SIZE} rows to be processed.")
        df = df.head(PILOT_SIZE)

    if BATCH_LIMIT:
        print(f"Batch Limit Active: {BATCH_LIMIT} new rows to be processed.")

    run_foldx_process(df)