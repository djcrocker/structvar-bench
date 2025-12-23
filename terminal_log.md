# Checking `verify_integrity.py`
```
dcrocker@Daniels-MacBook-Pro-684 structvar-bench % python3 src/data/verify_integrity.py
Starting Integrity Verification

Checking ClinVar Data: /Users/dcrocker/Projects/structvar-bench/data/raw/variant_summary.txt
/Users/dcrocker/Projects/structvar-bench/src/data/verify_integrity.py:20: DtypeWarning: Columns (18) have mixed types. Specify dtype option on import or set low_memory=False.
  df = pd.read_csv(CLINVAR_PATH, sep='\t') # Tab-separated
SUCCESS: File loaded successfully. Total Rows: 8367041
SUCCESS: All required columns found: ['Name', 'Chromosome', 'ClinicalSignificance', 'ReviewStatus']

Snapshot of first 5 variants:
                                                Name  ...                                       ReviewStatus
0  NM_014855.3(AP5Z1):c.80_83delinsTGCTGTAAACTGTA...  ...  criteria provided, multiple submitters, no con...
1  NM_014855.3(AP5Z1):c.80_83delinsTGCTGTAAACTGTA...  ...  criteria provided, multiple submitters, no con...
2     NM_014855.3(AP5Z1):c.1413_1426del (p.Leu473fs)  ...                     no assertion criteria provided
3     NM_014855.3(AP5Z1):c.1413_1426del (p.Leu473fs)  ...                     no assertion criteria provided
4       NM_014630.3(ZNF592):c.3136G>A (p.Gly1046Arg)  ...                     no assertion criteria provided

[5 rows x 4 columns]
Checking AlphaFold Structure with test file: /Users/dcrocker/Projects/structvar-bench/data/raw/alphaFold_human/AF-A0A0A0MRZ7-F1-model_v6.cif.gz
SUCCESS: CIF.gz file parsed successfully.
SUCCESS: B-factor (pLDDT) found. Value for first atom: 40.59
VALIDATION: pLDDT score is within expected range (0-100).

Verification Complete
```

# Checking `build_cohort.py`
```
dcrocker@Daniels-MacBook-Pro-684 structvar-bench % python3 src/data/build_cohort.py

04:37:55: Loading ClinVar data...
04:40:29: ClinVar data loaded.

Initial ClinVar rows: 8367041
Rows after Assembly 'GRCh38' filter (1): 4151373
Rows after SNV filter (2): 3843729
Rows after ReviewStatus filter (3): 3754105
Rows after Missense Extraction (4A): 2321577
Rows after Nonsense & Silent Eliminations (4B): 2240579
Rows after Class filter (Pathogenic/Benign): 186806

04:42:39: Loading UniProt Mapping file...
Rows after Mapping to UniProt IDs (5): 187117

04:42:39: Saving mapped cohort to /Users/dcrocker/Projects/structvar-bench/data/processed/cohort_mapped.csv
Saved final cohort.
```

# Checking `filter_structures.py`
```
dcrocker@Daniels-MacBook-Pro-684 structvar-bench % python3 src/data/filter_structures.py

04:46:15: Loading /Users/dcrocker/Projects/structvar-bench/data/processed/cohort_mapped.csv
04:46:15: Scanning /Users/dcrocker/Projects/structvar-bench/data/raw/alphaFold_human for structure files...
04:46:15: Found 47173 available structure files.
Found 20550 unique matching structure files.
Dropped 1101 rows missing structure files.

04:46:15: Checking pLDDT scores...
Processing 186000/186016...
15:46:42: Processing Complete.
Final Cohort Size: 99592

15:46:43: Saving final cohort to /Users/dcrocker/Projects/structvar-bench/data/processed/cohort_filtered.csv
```

# Checking `run_foldx.py` (showcasing resume feature)
```
dcrocker@Daniels-MacBook-Pro-684 structvar-bench % python3 src/local/run_foldx.py
--- FoldX Pipeline (Pilot Mode: True) ---
Sampling top 5 rows for testing...

19:05:32: Starting processing for 4 unique proteins...

19:05:32: [1/4] Processing Q30201 (1 mutations)...
  > Repaired structure found. Using cached version.

19:05:32: Protein repaired. Running mutations...
  > Mutating QA283P...
  > Done. ddG: 1.10636

19:05:38: [2/4] Processing Q86XE5 (2 mutations)...
  > Repaired structure found. Using cached version.

19:05:38: Protein repaired. Running mutations...
  > Mutating GA287V...
  > Done. ddG: 6.66249
  > Mutating RA97C...
  > Done. ddG: 1.72833

19:05:50: [3/4] Processing Q96CU9 (1 mutations)...
  > Repaired structure found. Using cached version.

19:05:50: Protein repaired. Running mutations...
  > Mutating NA430S...
  > Done. ddG: 3.47226

19:05:54: [4/4] Processing Q9P2L0 (1 mutations)...
  > Repaired structure found. Using cached version.

19:05:54: Protein repaired. Running mutations...
  > Mutating AA864T...
  > FoldX BuildModel failed for AA864T. Skipping.

Processing complete. Saving results...
Results saved to /Users/dcrocker/Projects/structvar-bench/data/processed/cohort_with_ddg.csv
Processed 5 variants.
```