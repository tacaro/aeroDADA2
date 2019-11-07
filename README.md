# aeroDADA2
This repository is a submodule. It contains asvs, sequences, and metadata for 16S amplicon data collected from C-20A science flights in the troposphere and lower stratosphere. The original sequences are publicly available on NASA GeneLab: https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-170/

### Files
- README.md: This file!
- aeRo.R: R script implementation of DADA2 pipeline
- aeRo.Rproj: corresponding R studio project
- asv_table.csv: ASV table with tabulated ASV counts corresponding to each sampling sequence
- sample_data.csv: Metadata for each sampling sequence.
  - Sample.Name: Name of sample
  - Material Type: Extracted DNA for all samples
  - Sample.Description: Qualitative comments on sampling, flight path, negative control, etc.
  - Organism: Metagenomic data for all samples
  - Flight.Name: C20A flight location
  - Sample.Category: Ground sampling or atmosphere sampling
  - Collection.Date: Date of sampling
  - Protocol.REF: Sample concentration for all
  - is.neg: Whether the sample represents a negative control
- seqtab_nochim.csv: tabulated sequences with chimeras removed
- tax_table.csv: ASVs labeled with taxa annotated by silva132 database.
- plots: subfolder containing plots generated from DADA2 pipeline.
  - fnFsquality.tiff: Quality score values of forward reads
  - fnRsquality.tiff: Quality score values of reverse reads
  - errorF.tiff: error rates, observed and predicted, for forward reads
  - errorR.tiff: error rates, observed and predicted, for reverse reads
