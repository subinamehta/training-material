The third and the last proteogenomics workflow is for identifying the novel peptides using BlastP and to localize the peptides to its genomic coordinates. Inputs from both workflow 1 and 2 will be used in this workflow.
The inputs for this workflow are:
● Tabular file – “Peptides for BlastP analysis”
● Tabular file – “PeptideShaker_PSM”
● Mz to sqlite
● Genomic mapping sqlite

All the files to run this workflow can be obtained from the second workflow output.Once the tabular output is created, we convert this tabular report into a FASTA file. This can be achieved by using the Tabular to FASTA convertion tool.

### Tabular to FASTA

> **Title column**: `1`
> **Sequence Column**:`2`


Now that we have the FASTA file, this is going to be subjected to BLAST-P analysis

### BLASTP (Basic Local Alignment Search Tool- proteins)

BLAST[https://blast.ncbi.nlm.nih.gov/Blast.cgi] is a web based tool used to compare biological sequences. BlastP, matches protein sequences against a protein database. More specifically, it looks at the amino acid sequence of proteins and can detect and evaluate the amount of differences between say, an experimentally derived sequence and all known amino acid sequences from a database. It can then find the most similar sequences and allow for identification of known proteins or for identification of potential peptides associated with novel proteoforms.

BlastP search is carried out with the PSM report (output from PeptideShaker). Before, BlastP analysis the “Peptides_for_Blast-P_analysis” is first converted from Tabular format to FASTA file format which can be easily read by the BlastP algorithm. This is done with the help of “Tabular to FASTA” conversion tool.
The short BlastP uses parameters for short peptide sequences (8-30 aas). Please use the rerun option to look at the parameters used.

>   **Protein query sequence(s)**: `Data input 'query' (fasta)`
>   **Subject database/sequences**: `Locally installed BLAST database`
>> - **Protein Blast database**: Select `nr_mouse_current`
>
>   **Type of BLAST**:`blastp-short - BLASTP optimized for queries shorter than 30 residues`
>
>   **Set expectation value cutoff**:`200000.0`
>
>   **Output format**: `Tabular(extended 25 columns)`
>
>   **Advanced Options**
>
>   **Filter out low complexity regions (with SEG)**: `No`
>
>   **Scoring matrix and gap costs**: `PAM30`
>
>   **Gap Costs**: `Existence: 9 Extension: 1`
>
>   **Maximum hits to show**: `1`
>
>   **Maximum number of HSPs (alignments) to keep for any single query-subject pair**:`1`
>
>   **Word size for wordfinder algorithm**: `2`
>
>   **Multiple hits window size**: `40`
>
>   **Minimum score to add a word to the BLAST lookup table**: `11`
>
>   **Composition-based statistics**: `0 (No composition)`
>
>   **Should the query and subject defline(s) be parsed?**: `No`
>
>   **Restrict search of database to a given set of ID's**:`No restriction`
>
>>       This feature provides a means to exclude ID's from a BLAST database search. 
>>       The expectation values in the BLAST results are based upon the sequences actually 
>>       searched, and not on the underlying database. Note this cannot be used when 
>>       comparing against a FASTA file.
>
>   **Minimum query coverage per hsp (percentage, 0 to 100)**: `0`
>
>   **Compute locally optimal Smith-Waterman alignments**:`No`

Once BlastP search is performed, it provides with a tabular output containing “Novel peptides”. Now this output is further processed by comparing the Novel Peptide output with the PSM report for selecting only distinct peptides which pass these criteria.
