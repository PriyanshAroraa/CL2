# CL2 - Bioinformatics Practicals

This repository contains bioinformatics practical assignments covering DNA sequence analysis, RNA-seq differential expression, protein structure prediction, molecular docking, machine learning for genomics, and GWAS analysis.

## Repository Structure

All files are organized in the `initial/` folder:

```
CL2/
├── initial/
│   ├── Code Files (.txt)
│   ├── Notes Files (.txt)
│   ├── Jupyter Notebooks (.ipynb)
│   └── Data Files (.txt, .csv, .fasta, .pdb, etc.)
└── README.md
```

## Contents

### Practical 1: DNA Sequence Analysis
- **Code**: `initial/BI1_DNA_Sequence_Analysis_Motif_GC_Content_ORF_Detection.txt`
- **Notes**: `initial/BI1_Notes_DNA_Sequence_Analysis.txt`
- **Notebook**: `initial/BI1.ipynb`
- **Description**: Analyze DNA sequences to find motifs, calculate GC content, and identify coding regions (ORFs)

### Practical 2: RNA-Seq Differential Expression Analysis
- **Code**: `initial/BI2_RNA_Seq_Differential_Expression_Analysis_DESeq2.txt`
- **Notes**: `initial/BI2_Notes_RNA_Seq_Differential_Expression.txt`
- **Notebook**: `initial/BI2.ipynb`
- **Description**: Perform differential gene expression analysis using DESeq2 on RNA-seq data

### Practical 3: Protein 3D Structure Prediction
- **Code**: `initial/BI3_Protein_3D_Structure_Prediction_Visualization.txt`
- **Notes**: `initial/BI3_Notes_Protein_3D_Structure.txt`
- **Notebook**: `initial/BI3.ipynb`
- **Description**: Read protein sequences from FASTA files and visualize 3D protein structures from PDB files

### Practical 4: Molecular Docking and Virtual Screening
- **Code**: `initial/BI4_Molecular_Docking_Virtual_Screening.txt`
- **Notes**: `initial/BI4_Notes_Molecular_Docking.txt`
- **Notebook**: `initial/BI4.ipynb`
- **Description**: Perform molecular docking simulations and virtual screening to identify potential drug candidates

### Practical 5: Machine Learning for Genomic Data Classification
- **Code**: `initial/BI5_Machine_Learning_Genomic_Data_Classification.txt`
- **Notes**: `initial/BI5_Notes_Machine_Learning_Genomics.txt`
- **Notebook**: `initial/BI5.ipynb`
- **Description**: Apply Random Forest and SVM algorithms to classify genomic sequences based on k-mer features

### Practical 6: GWAS Analysis for Crop Traits
- **Code**: `initial/BI6_GWAS_Crop_Traits_Genetic_Markers.txt`
- **Notes**: `initial/BI6_Notes_GWAS_Crop_Traits.txt`
- **Notebook**: `initial/BI6.ipynb`
- **Description**: Perform Genome-Wide Association Study (GWAS) to identify genetic markers associated with crop traits

## Data Files

All data files are located in the `initial/` folder:
- `human datset.txt` - Human DNA sequences dataset
- `Phenos.csv` - Phenotypic data for crop traits
- `sequence.fasta` - Protein sequence file
- `5RGU.pdb` - Protein structure file
- `pdb3rgk.ent` - Protein structure file
- `protein_ligand_docking_results.csv` - Docking results
- `results.sdf` - Structure data file

## Usage

Each practical has:
1. A `.txt` code file that can be copied into Python scripts
2. A notes file with comprehensive Q&A for viva preparation
3. A Jupyter notebook with the original implementation

Simply navigate to the `initial/` folder and copy the code from the `.txt` files into Python files and run them. Make sure to have the required data files in the same directory.

## Requirements

- Python 3.x
- Biopython
- pandas
- numpy
- matplotlib
- seaborn
- scikit-learn
- pydeseq2
- rdkit
- py3Dmol
- statsmodels

## Installation

```bash
pip install biopython pandas numpy matplotlib seaborn scikit-learn pydeseq2 rdkit py3dmol statsmodels
```

## License

This repository contains educational materials for bioinformatics practicals.

