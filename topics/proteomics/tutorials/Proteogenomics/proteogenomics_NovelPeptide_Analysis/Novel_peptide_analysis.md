The third and the last proteogenomics workflow is for identifying the novel peptides using BlastP and to localize the peptides to its genomic coordinates. Inputs from both workflow 1 and 2 will be used in this workflow.
The inputs for this workflow are:
> - Tabular file – “Peptides from BlastP analysis”
> - Tabular file – “PeptideShaker_PSM”
> - Mz to sqlite
> - Genomic mapping sqlite

All the files to run this workflow can be obtained from the second workflow output.Once the tabular output is created, we convert this tabular report into a FASTA file. This can be achieved by using the Tabular to FASTA convertion tool.


Once Blast-P search is performed, it provides with a tabular output containing “Novel peptides”. Now this output is further processed by comparing the Novel Peptide output with the PSM report for selecting only distinct peptides which pass these parameters:

### Query tabular (Extract novel peptides after BlastP)
1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - (a)**Database Table**: Click on `+ Insert Database Table`:
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `blast`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**:
>`qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,sallseqid,score,nident,positive,gaps,ppos,qframe,sframe,qseq,sseq,qlen,slen,salltitles`
>    - **Only load the columns you have named into database**: `Yes` 
>
>    Section **Table Index**:
>    - **Table Index**: `No`
>    - **Index on Columns**: `id`
>  
>     - (c) **Database Table**: Click on `+ Insert Database Table`:
>    Section **Filter Dataset Input**
>      - **Filter Tabular input lines**
>      - **Filter by**:  `skip leading lines`
>      - **Skip lines**: `1`
>
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**: `ID,Proteins,Sequence,AAs_Before,AAs_After,Position,Modified_Sequence,Variable_Modifications,Fixed_Modifications,Spectrum_File,Spectrum_Title,Spectrum_Scan_Number,RT,mz,Measured_Charge,Identification_Charge,Theoretical_Mass,Isotope_Number,Precursor_mz_Error_ppm,Localization_Confidence,Probabilistic_PTM_score,Dscore,Confidence,Validation`
>    - **Only load the columns you have named into database**: `Yes` 
>
>
>    - **Save the sqlite database in your history**: `No`
>
>        > ### {% icon tip %} Tip
>        >
>        > * **Query Tabular** can also use an existing SQLite database. Activating `Save the sqlite database in your history`
>        > will store the created database in the history, allowing to reuse it directly.
>        >
>        {: .tip}
>
>    - **SQL Query to generate tabular output**:
>
>>       SELECT distinct psm.*
>>       FROM psm join blast on psm.Sequence = blast.qseqid
>>       WHERE blast.pident < 100 OR blast.gapopen >= 1 OR blast.length < blast.qlen
>>       ORDER BY psm.Sequence, psm.ID 
       
>
>    - **include query result column headers**: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>

### Query tabular (Distinct Peptides from the list of Novel peptides)
1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - (a)**Database Table**: Click on `+ Insert Database Table`:
>
>    Section **Filter Dataset Input**
>    - **Filter Tabular input lines**
>      - Filter by:  `skip leading lines`
>      - Skip lines: `1`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**:`ID,Proteins,Sequence`
>    - **Only load the columns you have named into database**: `Yes` 
>
>    - **SQL Query to generate tabular output**:
>
>>       select distinct Sequence from psm
>
>    - **include query result column headers**: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>

### MVP

### Uploading only those peptides which have good lorikeet spectra

The next tool in the workflow is the Peptide genomic coordinate tool which takes the novel peptides as the input along with the mztosqlite file and the genomic mapping sqlite file (obtained during creation of the database). This tool helps create a bed file with the genomic coordinate information of the peptides based on the sqlite files. 

### Peptide genomic Coordinate


### Peppointer

### Query tabular( Final Summary)
1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - (a)**Database Table**: Click on `+ Insert Database Table`:
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `bed_pep_pointer`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**:`chrom,start,end,peptide,score,strand,annot`
>    - **Only load the columns you have named into database**: `No` 
>
>
>     - (b) **Database Table**: Click on `+ Insert Database Table`:
>    Section **Filter Dataset Input**
>      - **Filter Tabular input lines**
>      - **Filter by**:  `skip leading lines`
>      - **Skip lines**: `1`
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**: `ID,Proteins,Sequence,AAs_Before,AAs_After,Position,Modified_Sequence,Variable_Modifications,Fixed_Modifications,Spectrum_File,Spectrum_Title,Spectrum_Scan_Number,RT,mz,Measured_Charge,Identification_Charge,Theoretical_Mass,Isotope_Number,Precursor_mz_Error_ppm,Localization_Confidence,Probabilistic_PTM_score,Dscore,Confidence,Validation`
>    - **Only load the columns you have named into database**: `No` 
>
>    - **SQL Query to generate tabular output**:
>
>>       SELECT psm.Sequence as PeptideSequence, count(psm.Sequence) as SpectralCount, psm.Proteins as  Proteins,bed_pep_pointer.chrom as Chromosome, bed_pep_pointer.start as Start, bed_pep_pointer.end as End, bed_pep_pointer.strand as Strand, bed_pep_pointer.annot as Annotation, bed_pep_pointer.chrom||':'||bed_pep_pointer.start||'-'||bed_pep_pointer.end as GenomeCoordinate,'https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&position='||bed_pep_pointer.chrom||'%3A'||bed_pep_pointer.start||'-'||bed_pep_pointer.end as UCSC_Genome_Browser 
FROM psm 
INNER JOIN bed_pep_pointer on bed_pep_pointer.peptide = psm.Sequence 
GROUP BY psm.Sequence
>
>    - **include query result column headers**: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>
