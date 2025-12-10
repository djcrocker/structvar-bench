# StructVar-Bench: A Multi-Modal Dataset for Predicting Missense Variant Pathogenicity using AlphaFold Structures
Clinicians often find genetic mutations (variants) in patients but don’t know if they are harmless or disease-causing. AlphaFold gives us the structure, while tools like FoldX can calculate how a mutation changes stability. This project focuses on creating a clean, pre-calculated dataset of these structural energy changes to help train AI models.

# 📁 Resources Used and Downloads
- ClinVar is a public database of all human genetic variants.
    - For this project, we use the ClinVar txt summary, a clean, tab-delimited text file that already has headers outlined: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ -> variant_summary.txt.gz
    - For future work, the full ClinVar XML dump may be used for more thorough labeling: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/
- UnitprotKB is a public database of protein sequence and functional information.
    - For reviewed human entries, do the following:
        - Go to UniprotKB: https://www.uniprot.org/uniprotkb
        - Search for reviewed humans: `(taxonomy_id:9606) AND (reviewed:true)` -> https://www.uniprot.org/uniprotkb?query=*%28taxonomy_id%3A9606%29+AND+%28reviewed%3Atrue%29
        - Download FASTA (canonical) for sequences
        - Download TSV
            - Click customize columns and select PDB, AlphaFoldDB, and RefSeq for cross-referencing
    - For all UnitprotKB reviewed resources: https://www.uniprot.org/uniprotkb
- AlphaFold PDB provides access to hundreds of millions of protein structure predictions: https://alphafold.ebi.ac.uk/download
    - Here, we only care about the human proteins, so: https://alphafold.ebi.ac.uk/download -> UP000005640_9606_HUMAN_v4.tar

# 🗺️ File Structure
```
src
 ├── data/
      ├── raw/
           ├── variant_summary.txt # Ground truth labels
           ├── UP000005640_9606_HUMAN_v4/ # 3D structures in .pdb.gz format
           ├── human_reviewed.fasta # Sequences
           ├── human_id_mapping.tsv # Metadata table with RefSeq and AlphaFold columns         
```