# CL2 - Bioinformatics Practicals

This repository contains bioinformatics practical assignments covering DNA sequence analysis, RNA-seq differential expression, protein structure prediction, molecular docking, machine learning for genomics, and GWAS analysis.

## Contents

### Practical 1: DNA Sequence Analysis
- **Code**: `BI1_DNA_Sequence_Analysis_Motif_GC_Content_ORF_Detection.txt`
- **Notes**: `BI1_Notes_DNA_Sequence_Analysis.txt`
- **Description**: Analyze DNA sequences to find motifs, calculate GC content, and identify coding regions (ORFs)

### Practical 2: RNA-Seq Differential Expression Analysis
- **Code**: `BI2_RNA_Seq_Differential_Expression_Analysis_DESeq2.txt`
- **Notes**: `BI2_Notes_RNA_Seq_Differential_Expression.txt`
- **Description**: Perform differential gene expression analysis using DESeq2 on RNA-seq data

### Practical 3: Protein 3D Structure Prediction
- **Code**: `BI3_Protein_3D_Structure_Prediction_Visualization.txt`
- **Notes**: `BI3_Notes_Protein_3D_Structure.txt`
- **Description**: Read protein sequences from FASTA files and visualize 3D protein structures from PDB files

### Practical 4: Molecular Docking and Virtual Screening
- **Code**: `BI4_Molecular_Docking_Virtual_Screening.txt`
- **Notes**: `BI4_Notes_Molecular_Docking.txt`
- **Description**: Perform molecular docking simulations and virtual screening to identify potential drug candidates

### Practical 5: Machine Learning for Genomic Data Classification
- **Code**: `BI5_Machine_Learning_Genomic_Data_Classification.txt`
- **Notes**: `BI5_Notes_Machine_Learning_Genomics.txt`
- **Description**: Apply Random Forest and SVM algorithms to classify genomic sequences based on k-mer features

### Practical 6: GWAS Analysis for Crop Traits
- **Code**: `BI6_GWAS_Crop_Traits_Genetic_Markers.txt`
- **Notes**: `BI6_Notes_GWAS_Crop_Traits.txt`
- **Description**: Perform Genome-Wide Association Study (GWAS) to identify genetic markers associated with crop traits

## Data Files

- `human datset.txt` - Human DNA sequences dataset
- `Phenos.csv` - Phenotypic data for crop traits
- `sequence.fasta` - Protein sequence file
- `5RGU.pdb` - Protein structure file
- `pdb3rgk.ent` - Protein structure file
- `protein_ligand_docking_results.csv` - Docking results
- `results.sdf` - Structure data file

## Jupyter Notebooks

- `BI1.ipynb` - Practical 1 notebook
- `BI2.ipynb` - Practical 2 notebook
- `BI3.ipynb` - Practical 3 notebook
- `BI4.ipynb` - Practical 4 notebook
- `BI5.ipynb` - Practical 5 notebook
- `BI6.ipynb` - Practical 6 notebook

## Usage

Each practical has:
1. A `.txt` code file that can be copied into Python scripts
2. A notes file with comprehensive Q&A for viva preparation

Simply copy the code from the `.txt` files into Python files and run them. Make sure to have the required data files in the same directory.

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

