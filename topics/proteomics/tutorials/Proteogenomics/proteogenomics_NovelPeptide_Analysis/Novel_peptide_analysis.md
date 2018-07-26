The third and the last proteogenomics workflow is for identifying the novel peptides using BlastP and to localize the peptides to its genomic coordinates. Inputs from both workflow 1 and 2 will be used in this workflow.
The inputs for this workflow are:
> - Tabular file – “Peptides for BlastP analysis”
> - Tabular file – “PeptideShaker_PSM”
> - Mz to sqlite
> - Genomic mapping sqlite

All the files to run this workflow can be obtained from the second workflow output.Once the tabular output is created, we convert this tabular report into a FASTA file. This can be achieved by using the Tabular to FASTA convertion tool.


Once Blast-P search is performed, it provides with a tabular output containing “Novel peptides”. Now this output is further processed by comparing the Novel Peptide output with the PSM report for selecting only distinct peptides which pass these parameters:

### Query tabular (Extract novel peptides after BlastP)

### Query tabular (Distinct Peptides from the list of Novel peptides)

### MVP

### Uploading only those peptides which have good lorikeet spectra

### Peptide genomic Coordinate

### Peppointer

### Query tabular( Final Summary)
