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
>        FROM psm join blast on psm.Sequence = blast.qseqid
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
>    - **Tabular Dataset for Table**: Uniprot
>
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**:`prot`
>
> _**Table Index**_:
>    -**Table Index**: `No`
>    -**Index on Columns**: `Prot`
>
>    - (b) **Database Table**: Click on `+ Insert Database Table`:
>    Section **Filter Dataset Input**
>    - **Filter Tabular input lines**
>      - Filter by:  `skip leading lines`
>      - Skip lines: `1`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psms`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**:`id,Proteins,Sequence`
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
>    Add another filter tabular input lines
>      - **Filter Tabular input lines**
>      - **Filter by**:  `select columns`
>      - **Enter column numbers to keep**: `1,2`
>
>   Add another filter tabular input lines
>      - **Filter Tabular input lines**
>      - **Filter by**:  `normalize list columns,replicate rows for each item in the list`
>      - **Enter column numbers to normalize**: `2`
>      - **List item delimiter in column**: `,`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `prots`    
>    - **Use first line as column names** : `No`
>    - **Specify Column Names (comma-separated list)**:`id,prot`
>    - **Only load the columns you have named into database**: `Yes` 
>
> _**Table Index**_:
>    -**Table Index**: `No`
>    -**Index on Columns**: `prot, id`
>        > ### {% icon comment %} Comment
>        >
>        > By default, table columns will be named: c1,c2,c3,...,cn (column names for a table must be unique).
>        > You can override the default names by entering a comma separated list of names, e.g. `,name1,,,name2`
>        > would rename the second and fifth columns.
>        >
>        > Check your input file to find the settings which best fits your needs.
>        {: .comment}
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
>>       SELECT psms.* 
>>       FROM psms 
>>       WHERE psms.id NOT IN 
>>       (SELECT distinct prots.id 
>>       FROM prots JOIN uniprot ON prots.prot = uniprot.prot) 
>>       ORDER BY psms.id
>
>
>    - **include query result column headers**: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>
