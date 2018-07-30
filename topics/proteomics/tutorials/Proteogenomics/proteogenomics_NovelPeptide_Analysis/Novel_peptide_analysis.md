The third and the last proteogenomics workflow is for identifying the novel peptides using BlastP and to localize the peptides to its genomic coordinates. Inputs from both workflow 1 and 2 will be used in this workflow.
The inputs for this workflow are:
> - Tabular file – “Peptides from BlastP analysis”
> - Tabular file – “PeptideShaker_PSM”
> - Mz to sqlite
> - Genomic mapping sqlite

All the files to run this workflow can be obtained from the second workflow output.Once the tabular output is created, we convert this tabular report into a FASTA file. This can be achieved by using the Tabular to FASTA convertion tool.


Once Blast-P search is performed, it provides with a tabular output containing “Novel peptides”. Now this output is further processed by comparing the Novel Peptide output with the PSM report for selecting only distinct peptides which pass these parameters:

### Query tabular (Extract novel peptides after BlastP)

### Query tabular (Distinct Peptides from the list of Novel peptides)

### MVP

### Uploading only those peptides which have good lorikeet spectra

The next tool in the workflow is the Peptide genomic coordinate tool which takes the novel peptides as the input along with the mztosqlite file and the genomic mapping sqlite file (obtained during creation of the database). This tool helps create a bed file with the genomic coordinate information of the peptides based on the sqlite files. 

### Peptide genomic Coordinate


### Peppointer

### Query tabular( Final Summary)
