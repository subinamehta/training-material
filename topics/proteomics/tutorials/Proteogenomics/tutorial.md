---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: Proteogenomics_RNAseq_db_creation
---

# Introduction

**Proteogenomics** is a combination of proteomics, genomics and transcriptomics data to identify peptides and to understand protein-level evidence of gene expression. In this tutorial, we will create a protein database (FASTA) using RNA-sequencing files (FASTQ) and then perform database searching of the resulting FASTA file with the MS/MS data to identify novel peptides. Then, we will assign the genomic coordinates and annotations for these novel peptides as well as perform visualization of the data. 


<img src="training-material/topics/proteomics/images/potential_novel_publication.png" width=100%>

Proteogenomics most commonly integrates **RNA-Seq** data for generating customized protein sequence databases with mass spectrometry-based proteomics data, which are matched to these databases to identify novel protein sequence variants. (Cancer Res. (2017); 77(21):e43-e46. doi: <a target="_blank" href="https://doi.org/10.1158/0008-5472.CAN-17-0331">10.1158/0008-5472.CAN-17-0331</a>).


<img src="../../../images/workflow_objective1.png" width=100%>


## Part I


In this tutorial, the proteins and the total RNA were obtained from the early development of B-cells from mice. It was obtained at two developmental stages of B-cells: *Ebf1* -/- pre-pro-B and *Rag2* -/- pro-B. Please refer to the original study for details: [Heydarian, M. et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4276347/).


### Agenda
>
> In this tutorial, we will deal with:
>

> - _Pretreatments / Data upload_ 

> - _Alignment of RNA-sequencing data with reference genome_

> - _Creating FASTA files with SAVs, indels and splice junctions_

> - _Merging databases_


# Pretreatments

## Data upload

There are many ways to upload your data. Three among these are:

*   Uploading the files from your computer
*   Using a direct link
*   Importing from the data library if your instance provides the files

In this tutorial, we will get the data from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1302055.svg)](https://doi.org/10.5281/zenodo.1302055)

### Hands-on data upload and organization(:thumbsup: hands-on training)
>
> 1. Create a new history and name it something meaningful (e.g. *Proteogenomics DB creation*)
> 2. Import the two FASTQ files and the GTF file from Zenodo.

>>       Tip: Importing data via links:
>>          * Copy the link location 
>>          * Open the Galaxy Upload Manager
>>          * Select "Paste/Fetch Data"
>>          * Paste the link into the text field. You can add multiple links, each on a separate line.
>>          * Press "Start". By default, Galaxy uses the link as the dataset name.
>
>>       **Comments**: Rename the datasets with more descriptive names
>    
>
# Analysis

The first workflow focuses on creating a **FASTA** database created from RNA-seq data. There are two outputs from this workflow: (1) a **sequence database** consisting of variants and known reference sequences and (2) mapping files containing **genomic** and **variant** mapping data.


<img src="../../../images/database_creation.png" width=100%>

The first part of the workflow deals with creating the FASTA file containing sequences with single amino acid variants (SAVs), insertions and deletions (indels). The second part of the workflow helps in creating a FASTA file with transcript assemblies (splicing variants).

### Aligning FASTQ files on the human genome

The first tool in the workflow is the [**HISAT2**](http://ccb.jhu.edu/software/hisat) alignment tool. It maps next-generation sequence (NGS) reads to the reference genome. To run this tool, two input files are required: an RNA-seq file (.FASTQ) and a reference genome file in Gene transfer format (GTF). This .gtf file is obtained from the Ensembl database. When successful, this tool outputs a .bam file (binary version of a SAM: **S**equence **A**lignment/**M**ap).

#### HISAT2

> ### Hands-on: HISAT2
>
>  1. **HISAT2** :

> Parameters:

>    - **Source for the reference genome**: `Use a built-in genome` mm10
>    - **Single-end or paired-end reads**: `Single end` 
>    - **Input FASTQ files**: `FASTQ_ProB.fastqsanger`
>    - **Specify strand information**: `Unstranded`

>>     Tip: For paired inputs, select the paired end reads.

>   - **Section** **Summary Options**:
>    `Select default parameters`
>
>   - **Section** **Advanced Options**:
>    `Select default parameters`
    
>>      **Comments**: 
>>       Note that if your reads are from a stranded library, you need to choose the appropriate 
>>       setting for "Specify strand information" mentioned above. For single-end reads, use 'F' or 'R'. 
>>       'F' means that a read corresponds to a transcript. 'R' means that a read corresponds to the reverse 
>>       complemented counterpart of a transcript. For paired-end reads, use either 'FR' or 'RF'. With this 
>>       option being used, every read alignment will have an XS attribute tag: '+' means a read 
>>       belongs to a transcript on the positive '+' strand of the genome. '-' means a read belongs to a  
>>       transcript on the negative '-' strand of the genome. (TopHat has a similar option, --library--type option, 
>>       where fr - first strand corresponds to R and RF; fr - second strand corresponds to F and FR.)
 
> Once all parameters are selected, click **Execute**.



#### FreeBayes

[FreeBayes]( https://github.com/ekg/freebayes) is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.


<img src="../../../images/variant_calling.png" width=100%>


>>     **Comments**: Provided some BAM dataset(s) and a reference sequence, FreeBayes will produce 
>>     a VCF dataset describing SNPs, indels, and complex variants in samples in the input 
>>     alignments. By default, FreeBayes will consider variants supported by at least two observations  
>>     in a single sample (-C) and also by at least 20% of the reads from a single sample (-F). These 
>>     settings are suitable for low to high depth sequencing in haploid and diploid samples, but   
>>     users working with polyploid or pooled samples may adjust them depending on the characteristics 
>>     of their sequencing data.
>>
>>     FreeBayes is capable of calling variant haplotypes shorter than a read length where multiple 
>>     polymorphisms segregate on the same read. The maximum distance between polymorphisms phased in
>>     this way is determined by the --max-complex-gap, which defaults to 3bp. In practice, this can 
>>     comfortably be set to half the read length. Ploidy may be set to any level (-p), but by default 
>>     all samples are assumed to be diploid. FreeBayes can model per-sample and per-region variation  
>>     in copy-number (-A) using a copy-number variation map.
>>
>>     FreeBayes can act as a frequency-based pooled caller and can describe variants and haplotypes in 
>>     terms of observation frequency rather than called genotypes. To do so, use --pooled-continuous  
>>     and set input filters to a suitable level. Allele observation counts will be described by AO   
>>     and RO fields in the VCF output.


> 1. **FreeBayes** :
>   - **Choose the source for the reference genome**: `Locally cached file`
>       - **Run in batch mode?**: `Run Individually`
>   - **BAM dataset**: `HISAT_Output.BAM`
>   - **Using reference genome**: `mm10`
>   - **Limit variant calling to a set of regions?**: `Do not Limit`
>   - **Choose parameter selection level**: `Simple diploid calling`
>   
>
>>    **Comments**: Galaxy allows five levels of control over FreeBayes options, provided by the Choose 
>>    parameter selection level menu option. These are: 
>>    1. Simple diploid calling: The simplest possible FreeBayes application. Equivalent to using  
>>    FreeBayes with only a BAM input and no other parameter options.
>>
>>    2. Simple diploid calling with filtering and coverage: Same as #1 plus two additional options: 
>>     (1) -0 (standard filters: --min-mapping-quality 30 --min-base-quality 20 --min-supporting-allele 
>>     -qsum 0 --genotype-variant-threshold 0) and (2) --min-coverage.
>>
>>    3. Frequency-based pooled calling: This is equivalent to using FreeBayes with the following  
>>    options:--haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled
>>    -continuous  --report- monomorphic. This is the best choice for calling variants in mixtures   
>>    such as viral, bacterial, or organellar genomes.
>>
>>    4. Frequency-based pooled calling with filtering and coverage: Same as #3 but adds -0 and  
>>    --min-coverage like in #2.

> Complete list of all options: Gives you full control by exposing all FreeBayes options as Galaxy parameters.
>    
>
> Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=40>
>
>
>    
>   


#### CustomProDB

[CustomProDB]( http://dx.doi.org/10.1093/bioinformatics/btt543) generates custom protein FASTAs from exosome or transcriptome data. Once Freebayes creates the .vcf file, CustomProDB uses this file to generate a custom protein FASTA file from the transcriptome data. For this tool, we use Ensemble 89 mmusculus (GRm38.p5) (dbsnp142) as the genome annotation. We create three FASTA files from CustomProDB: (1) a variant FASTA file for short indels, (2) a Single Amino acid Variant (SAV) FASTA file, an SQLite database file (genome mapping and variant mapping) for mapping proteins to a genome and (3) an RData file for variant protein coding sequences.
The reference protein set can be filtered by transcript expression level (RPKM calculated from a BAM file), and variant protein forms can be predicted based on variant calls (SNPs and indels reported in a VCF file).

>>      **Comments**: Annotations CustomProDB depends on a set of annotation files (in RData format) to 
>>      create reference and variant protein sequences. Galaxy administrators can use the CustomProDB   
>>      data manager to create these annotations to make them available for users.
<img src="../../../images/CustomProDB.png" width=100%>

### Hands-on: Using CustomProDB to generate protein FASTAs from exosome or transcriptome data
>
> 1. **CustomProDB**:
>   - **Will you select a genome annotation from your history or use a built-in annotation?**: `Use built in genome annotation`
>   - **Using reference genome**: `Ensemble 89 mmusculus (GRm38.p5) (dbsnp142)`
>   - **BAM file**: `HISAT_Output.BAM`
>   - **VCF file**: `Freebayes.vcf`
>   - **Annotate SNPs with rsid from dbSNP**: `No`
>   - **Annotate somatic SNPs from COSMIC (human only)**: `No`
>   - **Transcript Expression Cutoff (RPKM)**: `1`
>   - **Create a variant FASTA for short insertions and deletions**: `Yes`
>   - **Create SQLite files for mapping proteins to genome and summarizing variant proteins**: `Yes`
>   - **Create RData file of variant protein coding sequences**: `Yes`


>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=30>
>
>>       **Comments**: Three FASTA files are created through the CustomProDB tool: a variant FASTA  
>>      file for short indels, a Single Amino acid Variant (SAV) FASTA file, an Sqlite file (genome   
>>      mapping and variant mapping) for mapping proteins to genome and an RData file for variant protein 
>>      coding sequences. Similar to the genomic mapping, a variant mapping file is also created from  
>>      CustomProDB. This SQLite file is also converted to tabular format and made SearchGUI-compatible. This 
>>      variant annotation file will be used to visualize the variants in the Multi-omics Visualization 
>>      Platform (in-house visualization platform developed by Galaxy-P senior developers).


#### StringTie

[StringTie](http://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional *de novo* assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. 

Its input can include not only the alignments of raw reads used by other transcript assemblers, but also alignments of longer sequences that have been assembled from those reads. To identify differentially expressed genes between experiments, StringTie's output can be processed by specialized software like Ballgown, Cuffdiff or other programs (DESeq2, edgeR, etc.).


1. **StringTie** 
> **StringTie transcript assembly and quantification**
>   - **Input mapped reads**: `FASTQ_ProB.BAM`
>   - **Specify strand information**: `Unstranded`
>   - **Use a reference file to guide assembly?**: `Use Reference GTF/GFF3`
>   - **Reference file**: `Use file from History`
>       - **GTF/GFF3 dataset to guide assembly**: `Mus_musculus.GRCm38.86.gtf`
>   - **Use Reference transcripts only?**: `No`
>   - **Output files for differential expression?**: `No additional output`
>   - **Output coverage file?**: `No`
>   - **Advanced Options**: `Default Parameters`


>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=30>
>    

>>     **Comments**:
>>     StringTie accepts a BAM (or SAM) file of paired-end RNA-seq reads, which must be  
>>     sorted by genomic location (coordinate position). This file contains spliced read alignments 
>>     and can be produced directly by programs such as HISAT2. We recommend using HISAT2 as it is a  
>>     fast and accurate alignment program. Every spliced read alignment (i.e. an alignment across  
>>     at least one junction) in the input BAM file must contain the tag XS to indicate the genomic  
>>     strand that produced the RNA from which the read was sequenced. Alignments produced by HISAT2  
>>     (when run with the --dta option) already include this tag, but if you use a different read mapper 
>>     you should check that this XS tag is included for spliced alignments.

>>     NOTE: be sure to run HISAT2 with the --dta option for alignment (under 'Spliced alignment options'),  
>>     or your results will suffer.
>>
>>     Also note that if your reads are from a stranded library, you need to choose the appropriate  
>>     setting for 'Specify strand information' above. As, if Forward (FR) is selected, StringTie will  
>>     assume the reads are from a --fr library, while if Reverse (RF) is selected, StringTie will 
>>     assume the reads are from a --rf library, otherwise it is assumed that the reads are from an 
>>     unstranded library (the widely-used, although now deprecated, TopHat had a similar --library
>>     -type option, where fr-firststrand corresponded to RF and fr-secondstrand corresponded to FR). 
>>     If you don't know whether your reads are from a stranded library or not, you could use the 
>>     tool 'RSeQC Infer Experiment' to try to determine.
>>     
>>     As an option, a reference annotation file in GTF/GFF3 format can be provided to StringTie. In  
>>     this case, StringTie will prefer to use these "known" genes from the annotation file, and for  
>>     the ones that are expressed it will compute coverage, TPM and FPKM values. It will also produce  
>>     additional transcripts to account for RNA-seq data that aren't covered by (or explained by) the 
>>     annotation. Note that if option -e is not used, the reference transcripts need to be fully covered 
>>     by reads in order to be included in StringTie's output. In that case, other transcripts assembled 
>>     from the data by StringTie and not present in the reference file will be printed as well.
>>     
>>
>>     We highly recommend that you provide annotation if you are analyzing a genome that is well
>>     annotated, such as human, mouse, or other model organisms.

#### GffCompare

[GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) compares and evaluates the accuracy of RNA-Seq transcript assemblers (Cufflinks, Stringtie). * collapse (merge) duplicate transcripts from multiple GTF/GFF3 files (e.g. resulted from assembly of different samples) * classify transcripts from one or multiple GTF/GFF3 files as they relate to reference transcripts provided in a annotation file (also in GTF/GFF3 format)

The original form of this program is also distributed as part of the Cufflinks suite, under the name ["CuffCompare"] (http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/). Most of the options and parameters of CuffCompare are supported by GffCompare, while new features will likely be added to GffCompare in the future.

1. **Using GffCompare to check assembled transcripts against a reference annotation**:

>   - **GTF inputs for comparison**`Stringtie_outut.gtf`
>   - **Use Reference Annotation**: `Mus_musculus.GRCm38.86.gtf`
>   - **Reference Annotation**: `Unstranded`

>   - **Ignore reference transcripts that are not overlapped by any input transfrags**: `No`
>   - **Ignore input transcripts that are not overlapped by any reference transcripts**: `No`
>   - **Use Sequence Data**: `No`
>   - **discard (ignore) single-exon transcripts**: `No`
>   - **Max. Distance for assessing exon accuracy**: `100`
>   - **Max distance for transcript grouping**: `100`
>   - **discard intron-redundant transfrags sharing 5'**: `No`

>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=30>
>     

>>      **Comments**:
>>      A notable difference between GffCompare and CuffCompare is that when a single query GTF/GFF file is 
>>      given as input along with a reference annotation (-r option), GFFCompare switches into "annotation mode" 
>>      and it generates a .annotated.gtf file instead of the .merged.gtf produced by CuffCompare with the same 
>>      parameters. This file has the same general format as CuffCompare's .merged.gtf file (with "class  
>>      codes" assigned to transcripts as per their relationship with the matching/overlapping reference, 
>>      transcript) but the original transcript IDs are preserved, so GffCompare can thus be used as a  
>>      simple way of annotating a set of transcripts.
>>
>>      Another important difference is that the input transcripts are no longer discarded when they  
>>      are found to be "intron redundant", i.e. contained within other, longer isoforms. CuffCompare  
>>      had the -G option to prevent collapsing of such intron redundant isoforms into their longer 
>>      "containers", but GffCompare has made this the default mode of operation (hence the -G option is
>>      no longer needed and is simply ignored when given).
>>      

#### Convert a GffCompare-annotated GTF file to BED for StringTie results

Convert a GffCompare annotated GTF file to BED format.

1. **GffCompare: compare assembled transcripts to a reference annotation**: `output from gff compare`
>
>  **GTF annotated by gffCompare**: 
>
>  **filter GffCompare class_codes to convert** 
>
>  - `j : Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript`
>  - `e : Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible   pre-mRNA fragment.`
>  - `i : A transfrag falling entirely within a reference intron`
>  - `p : Possible polymerase run-on fragment (within 2Kbases of a reference transcript)`
>  - `u : Unknown, intergenic transcript`
>
>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=30>
>     

#### Translate BED transcripts cDNA in 3frames or CDS
Translate transcripts from the input BED file into protein sequences.
 
> The genomic sequence:
1. may be supplied in an extra column in the BED input file
2. retrieved from a .2bit genomic reference file
3. retrieved from the Ensembl REST API for Ensembl transcripts

1. **Translate BED transcripts cDNA in 3frames or CDS** {% icon tool %}:
**A BED file with 12 columns**: `Convert GffCompare-annotated GTF to BED`
**Source for Genomic Sequence Data** `Locally cached File`
**Select reference 2bit file** `mm10`
**BED Filtering Options** `default`
**Translation Options ** `default`
**FASTA ID Options** `default`

>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=30>
>   

>>       **Comments**:
>>
>>       **INPUTS**
>>          BED file with at least the standard 12 columns
>>          Genome reference in .2bit format (optional)
>>
>>       **OUTPUTS**
>>
>>          FASTA of transcript translations
>>          BED with the genomic location of the translated protein. The added 13th column contains the 
>>          protein sequence.
>>
>>       **OPTIONS**
>>
>>       Feature translation
>>          - cDNA - three-frame translations of the cDNA sequences with an output for each sequence  
>>              between STOP codons
>>          - CDS - three-frame translations of CDS (coding sequence defined by thickStart and thickEnd  
>>              in the BED file)
>>       Translation filtering
>>          - can be trimmed to a Methionine start codon
>>          - can be split into peptides by an enzyme digestion
>>          - must exceed specified minimum length
>>       BED filtering
>>          - genomic regions
>>          - Ensembl biotype if the BED contains the 20 columns as retrieved from the Ensembl REST API
    
#### BED to protein map genomic location of proteins for MVP

Convert a BED format file of the proteins from a proteomics search database into a tabular format for the Multiomics Visualization Platform (MVP).

>   1. **A BED file with 12 columns, thickStart and thickEnd define protein coding region**: `Translate cDNA_minus_CDS`

>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon:<img src="../../../images/view_icon.png" height=30>
>   
>>       **Comments**:
>>      - The tabular output can be converted to a sqlite database using the Query_Tabular tool.
>>      - The sqlite table should be named: feature_cds_map.
>>      - The names for the columns should be: name, chrom, start, end, strand, cds_start, cds_end
>>      - This SQL query will return the genomic location for a peptide sequence in a protein    
>>        (multiply the amino acid position by 3 for the cds location)


### Creating FASTA Database:

The Protein Database Downloader tool is used to download the FASTA database from UniProt and cRAP database containing known/reference mouse proteins.

#### FASTA Merge Files and Filter Unique Sequences Concatenate FASTA database files together

- Concatenate FASTA database files together.
-  If the uniqueness criterion is "Accession and Sequence", only the first appearence of each unique sequence will appear in the output. Otherwise, duplicate sequences are allowed, but only the first appearance of each accession will appear in the output.
- The default accession parser will treat everything in the header before the first space as the accession.

> 1. **Run in batch mode?**: `Merge individual FASTAs (output collection if input is collection)`

> *Input FASTA File(s)* : ` Input Custom ProDB Fasta File output`
                        - ` 1.HISAT_Output.rpkm`
                        - ` 2.HISAT_Output.snv`
                        - ` 3.HISAT_Output.indel`

>  **How are sequences judged to be unique?**:`Accession and Sequence`
>
>  **Accession Parsing Regular Expression**: `^>([^ |]+).*$`

>   2. Click **Execute** and inspect the resulting files after they turn green with the **View data** icon: <img src="../../../images/view_icon.png" height=30>
>   

>>     **Comments**:
>>     The Regex Text Manipulation tool is used to manipulate the FASTA file to make it 
>>     SearchGUI-compatible. The “FASTA Merge Files and Filter Unique Sequences"  
>>     tool is used to merge the databases obtained from CustomProDB and "Translate BED" tool
>>     along with the UniProt and cRAP databases.
>


<img src="../../../images/Fasta_sequence.png" width=100%>

For visualization purposes we also use the concatenate tool to concatenate the genomic mapping with the protein mapping dataset. This output will be used for visualization in MVP to view the genomic coordinates of the variant peptide.


An SQLite database containing the genomic mapping SQLite, variant annotation and information from the protein mapping file is concatenated to form a single genomic mapping SQLite database later used as an input for the 'Peptide Genomic Coordinate' tool. For that we need to follow the steps below:

### SQLite to tabular for SQL query (for genomic mapping)

> The input for this tool is an existing SQLite database (genomic_mapping.sqlite from CustomProDB) and the outputs are the results of an SQL query exported to the history as a tabular file.


>  **Query**:
>
>      `SELECT pro_name, chr_name, cds_chr_start - 1, cds_chr_end,strand,cds_start - 1, cds_end
>      FROM genomic_mapping
>      ORDER BY pro_name, cds_start, cds_end`

We will subject the output to text manipulation so that the results are compatible with the Multiomics Visualization Platform.

### Column Regex Find And Replace (SearchGUI-compatible Protein Names Genomic Mapping)

> 2. Click **Execute**

This tool goes line by line through the specified input file and if the text in the selected column matches a specified regular expression pattern, it replaces the text with the specified replacement.

> **Select cells from**: `genomic_mapping_sqlite' (tabular)`
> **Using:** `column 1`
> Select insert check

>   **Check 1**: 
>>  **Find Regex:**
>>     `^(ENS[^_]+_\d+:)([ACGTacgt]+)>([ACGTacgt]+)\s*`
>>  **Replacement:**  
>>     `\1\2_\3`
>>  **Check 2**:
>>  **Find Regex:**
>>     `,([A-Z]\d+[A-Z])\s*`
>>  **Replacement:**  
>>     `.\1`
>>  **Check 3**:
>>  **Find Regex:**
>>    ` ^(ENS[^ |]*)\s*`
>>  **Replacement:**  
>>     \1

>   2. Click **Execute**
Once this step is complete we will concatenate the output from this tool with the "Bed to protein map" output.

### Concatenate multiple datasets
> Select the output from the previous tool with "Bed2protein_SJ_SAV_INDEL" output.

> Output will be the "Genomic_Protein_map"

### Query Tabular using SQLite (For genomic mapping)

Loads tabular datasets into an SQLite database.

1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: `Genomic_Protein_map`
>
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `feature_cds_map`
>    - **Specify Column Names (comma-separated list)**: `name,chrom,start,end,strand,cds_start,cds_end`
>
>    - **Only load the columns you have named into database**: `No`
>
>    Section **Table Index**:
>
>    - **This is a unique index**: `No`
>    - **Index on columns**: `name,cds_start,cds_end`
>
> Rename the output as **"genomic_mapping_sqlite"**
>
<img src="../../../images/genomic_mapping_file.png" width=100%>

> 2. Click **Execute** and inspect the query results file after it turns green:


### SQLite to tabular for SQL query (For variant annotations)

> The input for this tool is an existing SQLite database (genomic_mapping.sqlite from CustomProDB) and the outputs are the results of an SQL query that is exported to the history as a tabular file.


>  **Query**:
>
>      `SELECT var_pro_name,pro_name,cigar,annotation
>       FROM variant_annotation`

We will subject the output to text manipulation so that the results are compatible with the Multiomics Visualization Platform.

### Column Regex Find And Replace (SearchGUI compatible Protein Names Variant annotations)
This tool goes line by line through the specified input file and if the text in the selected column matches a specified regular expression pattern, it replaces the text with the corresponding specified replacement.

> **Select cells from**: `variant_annotations_sqlite' (tabular)`
> **Using:** `column 1`
> Select insert check

>   **Check 1**: 
>>  **Find Regex:**
>>     `^(ENS[^_]+_\d+:)([ACGTacgt]+)>([ACGTacgt]+)\s*`
>>  **Replacement:**  
>>     `\1\2_\3`
>>  **Check 2**:
>>  **Find Regex:**
>>     `,([A-Z]\d+[A-Z])\s*`
>>  **Replacement:**  
>>     `.\1`
>>  **Check 3**:
>>  **Find Regex:**
>>    ` ^(ENS[^ |]*)\s*`
>>  **Replacement:**  
>>     \1


### Query Tabular using SQLite SQL (For variant annotations)

>Loads tabular datasets into an SQLite database.

1. **Query Tabular** : Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: `variant_annotation`
>
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `variant_annotation`
>    - **Specify Column Names (comma-separated list)**: `name,reference,cigar,annotation`
>
>    - **Only load the columns you have named into database**: `No`
>
>    Section **Table Index**:
>
>    - **This is a unique index**: `No`
>    - **Index on columns**: `name,cigar`
>
> Rename the output as **"Variant_annotation_sqlitedb"**
>
>
> 2. Click **Execute** and inspect the query results file after it turns green:
>

>>      These SQLite databases, which contain the genomic mapping SQLite and variant 
>>      annotation information from the protein mapping file, will be utilized
>>      by MVP to visualize the genomic loci of any variant peptides.


<img src="../../../images/viewing_SNP_Variant_IGV.png" width=100%>

## Part II
---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: Proteogenomics_RNAseq_db_search
---

# Introduction

In this tutorial, we perform proteogenomic database searching using the Mass Spectrometry data. The inputs for performing the proteogenomic database searching are the MGF files and the FASTA database file. The FASTA database is obtained by running the first workflow “Uniprot_cRAP_SAV_indel_translatedbed.FASTA”. The second workflow focuses on performing database search of the peak list files (MGFs).
<img src="../../../images/second_workflow.png" width=100%>
> ### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> - _Pretreatments / Data upload_ 

> - _Database searching using SearchGUI and Peptide Shaker_

> - _Removing known peptides_

> - _Performing Blast-P analysis for obtaining novel proteoforms_

# Pretreatments


### Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (e.g. *Proteogenomics DB search*)
> 2. Import the three MGF MS/MS files and the FASTA sequence file from Zenodo.
>
>>       Tip: Importing data via links
>>          * Copy the link location 
>>          * Open the Galaxy Upload Manager
>>          * Select "Paste/Fetch Data"
>>          * Paste the link into the text field. You can add multiple links, each on a separate line.
>>          * Press Start. As default, Galaxy takes the link as name.
>>
>>       **Comments**:Rename the datasets to a more descriptive name
>    
>
> 3. Build a **Dataset list** for the three MGF files
>    - Click the **Operations on multiple datasets** check box at the top of the history panel
>       ![Operations on multiple datasets button]<img src="../../../images/dataset_list.png" width=10%>
>    - Check the three boxes next to the MGF files
>    - Click **For all selected...** and choose **Build dataset list**
>    - Ensure the three control samples are the only ones selected, and enter a name for the new collection (e.g. *MGF files*)
>    - Click **Create list** and exit by clicking again the dataset operations icon
>


# Analysis

## Match peptide sequences

The search database labelled `Uniprot_cRAP_SAV_indel_translatedbed.FASTA` is the input database that
will be used to match MS/MS to peptide sequences via a sequence database search. 

For this, the sequence database-searching program called [SearchGUI](https://compomics.github.io/projects/searchgui.html) will be used.The created dataset collection of the three *MGF files* in the history is used as the MS/MS input. We will walk through a number of these settings in order to utilize SearchGUI on these example MGF files.

#### SearchGUI

### Hands-on: SearchGUI
>
> 1. **SearchGUI**: Run **SearchGUI** with:
>    - **Protein Database**: `Uniprot_cRAP_SAV_indel_translatedbed.FASTA`(or however you named the `FASTA` file)
>    - **Input Peak lists (mgf)**: `MGF files` dataset collection.
>
>    >>     Tip: Select dataset collections as input
>    >
>    > * Click the **Dataset collection** icon on the left of the input field:
>    >
>    >      ![Dataset collection button]<img src="../../../images/dataset_list.png" width=10%>
>    > * Select the appropriate dataset collection from the list
>    {: .tip}
>
>    Section **Search Engine Options**:
>
>    - **B-Search Engines**: `X!Tandem`
>
>>      **Comments**:
>>    The section **Search Engine Options** contains a selection of sequence database searching
>>    programs that are available in SearchGUI. Any combination of these programs can be used for
>>    generating PSMs from MS/MS data. For the purpose of this tutorial, **X!Tandem** we will be 
>>    used.    
>
>    Section **Precursor Options**:
>   
>    **Enzyme**: `Trypsin`
>    **Maximum Missed Cleavages**: `2`
>    **Precursor Ion Tolerance Units**: `Parts per million (ppm)`
>    **Precursor Ion Tolerance**:` 10`
>    **Fragment Tolerance (Daltons)**: `0.05` (this is high resolution MS/MS data) 
>    **Minimum charge**:`2`
>    **Maximum charge**:`6`
>    **Forward Ion**: `b`
>    **Reverse Ion**:` y`
>    **Minimum Precursor Isotope** :`0`
>    **Maximum Precursor Isotope** :`1`
>
>    Section **Protein Modification Options**:
>
>    - **Fixed Modifications**: `Carbamidomethylation of C, ITRAQ-4Plex of K, ITRAQ-4Plex of Ntermini`
>    - **Variable modifications**: `Oxidation of M, ITRAQ-4Plex of Y`
>
>>       ### Tip: Search for options
>>       * For selection lists, typing the first few letters in the window will filter the
>>         available options.  
>
>    Section **Advanced Options**:
>    - **X!Tandem Options**: `Advanced`
>    - **X!Tandem: Quick Acetyl**: `No`
>    - **X!Tandem: Quick Pyrolidone**: `No`
>    - **X!Tandem: Protein stP Bias**: `No`
>    - **X!Tandem: Maximum Valid Expectation Value**: `100`
>
>    - leave everything else as default
>
> 2. Click **Execute**.
>


Once the database search is completed, the SearchGUI tool will output a file (called a
SearchGUI archive file) that will serve as an input for the next section, PeptideShaker.

>>    ### Comment
>>    Note that sequence databases used for metaproteomics are usually much larger than 
>>    the excerpt used in this tutorial. When using large databases, the peptide identification 
>>    step can take much more time for computation. In metaproteomics, choosing the optimal 
>>    database is a crucial step of your workflow, for further reading see 
>>    [Timmins-Schiffman et al (2017)](https://www.ncbi.nlm.nih.gov/pubmed/27824341).
>>



#### PeptideShaker

[PeptideShaker](https://compomics.github.io/projects/peptide-shaker.html) is a post-processing software tool that
processes data from the SearchGUI software tool. It serves to organize the Peptide-Spectral
Matches (PSMs) generated from SearchGUI processing and is contained in the SearchGUI archive.
It provides an assessment of confidence of the data, inferring proteins identified from the
matched peptide sequences and generates outputs that can be visualized by users to interpret
results. PeptideShaker has been wrapped in Galaxy to work in combination with SearchGUI
outputs.

>>     ### Comment
>>      There are a number of choices for different data files that can be generated using
>>      PeptideShaker. A compressed file can be made containing all information needed to 
>>      view them results in the standalone PeptideShaker viewer. A `mzidentML` file can 
>>      be created that contains all peptide sequence matching information and can be 
>>      utilized by compatible downstream software. Other outputs are focused on the inferred 
>>      proteins identified from the PSMs, as well as phosphorylation reports, relevant if 
>>      a phosphoproteomics experiment has been undertaken. 
>>


 ###  Hands-on: PeptideShaker
>
> 1. **PeptideShaker** : Run **PeptideShaker** with:
>   - **Compressed SearchGUI results**: The SearchGUI archive file
>   - **Specify Advanced PeptideShaker Processing Options**: `Default Processing Options`
>   - **Specify Advanced Filtering Options**: `Default Filtering Options`
>   - **Include the protein sequences in mzIdentML**: `No`
>   - **Output options**: Select the `PSM Report` (Peptide-Spectral Match) and the `Certificate of Analysis`
>
>>     ###  Comment
>>
>>     The **Certificate of Analysis** provides details on all the parameters
>>     used by both SearchGUI and PeptideShaker in the analysis. This can be downloaded from the
>>     Galaxy instance to your local computer in a text file if desired.
>      
>
> 2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button]<img src="../../../images/view_icon.png" width=15%>
>



A number of new items will appear in your history, each corresponding to the outputs selected
in the PeptideShaker parameters. Most relevant for this tutorial is the PSM report:

Scrolling at the bottom to the left will show the sequence for the PSM that matched to these
metapeptide entries. Column 3 is the sequence matched for each PSM entry. Every PSM is a
new row in the tabular output.

A number of new items will appear in your History, each corresponding to the outputs selected in the PeptideShaker parameters. The Peptide Shaker’s PSM report is used as an input for the BlastP analysis. Before performing BlastP analysis. The Query Tabular tool and few test manipulation tools are used to remove spectra that belongs to the reference proteins. The output tabular file “Peptides_for_Blast-P_analysis” will contain only those spectra that did not belong to any known proteins.

#### Creating SQLITE database using mz to sqlite

The mzidentml output from the Peptide shaker is converted into an sqlite database file by using the mz to sqlite tool. This sqlite output is used to open the Multi-omics visualization platform, wherein you can view the spectra of the peptides using Lorikeet parameters. To open the MVP viewer, click on the “Visualize in MVP Application” icon ( this will pop-open the interactive multi-omics viewer in a new window/tab)


### Hands-on: mz to sqlite: extract mzidentml ans associated proteomics datasets into a sqlite db

>
> 1. **mz to sqlite** : Run **mz to sqlite** with:
>
>    - **Proteomics identification files**: Click on `PeptideShaker_mzidentml`:
>    - **Proteomics Spectrum files**: `Mo_Tai_MGFs`
>    - **Proteomics Search Database Fasta**: `Uniprot_cRAP_SAV_indel_translatedbed.FASTA`
>
Click **Execute**
<img src="../../../images/mz2sqlite.png" width=50%>

The next step is to remove known peptides from the list of PSM's that we acquired from the Peptide shaker results. For that we need to perform some text manipulation steps to extract list of known peptides from the Uniprot and cRAP database.

### Merged Uniprot and cRAP database
The file named "Trimmed_ref_500_Uniprot_cRAP.fasta" is the trimmed version of Uniprot and cRAP database merged Fasta files.

This Fasta file will be subjected to few text manipulation steps in order to get the tabular file for the known peptides. The first step is to convert this FASTA file to tabular in order to proceed with text manipulation.

###  Hands-on: FASTA to Tabular
Convert these sequences:

> **Data input 'input' (fasta)**: `Trimmed_ref_500_Uniprot_cRAP.fasta`
> **How many columns to divide title string into?**: `2`
> **How many title characters to keep?**: `0`

The resultant tabular file will go through a series of text manipulation steps.


### Text Manipulation steps 

> 1. **Cut**
>    - **Cut Columns**: `c1`
>    - **Delimited by**: `Tab`
 
 Upon completion of this step you will have extracted C1 (column 1) from the input tabular file
 
> 2. **Convert**
>   - **Convert all**: `Whitespaces`
>   - **in Dataset** : `Data input 'input' (txt)`
>   - **Strip leading and trailing whitespaces**: `Yes`
>   - **Condense consecutive delimiters in one TAB**: `Yes`

This step will convert all the white spaces into different tabular column.

> 3. **Cut**
>    - **Cut Columns**: `c2`
>    - **Delimited by**: `Pipe`

This step will extract information in column2 separated by Pipe (|)

> 4. **Convert**
>   - **Convert all**: `Dots`
>   - **in Dataset** : `Data input 'input' (txt)`
>   - **Strip leading and trailing whitespaces**: `Yes`
>   - **Condense consecutive delimiters in one TAB**: `Yes`
This step will convert all the dots into different tabular column.

> 5. **Group**
>   - **Select data**: `input from above`
>   - **Group by column**: `1`


Now that we have the list of known peptides, the query tabular tool is used to move these reference pepides from the PSM report.


#### Query Tabular

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
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
>     ![Query Tabular output showing the peptides]
>
**The output from this step is that the resultant peptides would be those which doesn't belong in the Uniprot or cRAP database.The query tabular tool is used again to create a tabular output containing peptides ready for Blast P analysis.**
>
> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
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
>    - **Specify Column Names (comma-separated list)**:`id,Proteins,Sequence`
>    - **Only load the columns you have named into database**: `Yes` 
>
>    - **SQL Query to generate tabular output**:
>
>>       SELECT Sequence || ' PSM=' || group_concat(id,',') || ' length=' 
>>       || length(Sequence) as "ID",Sequence
>>       FROM  psm
>>       WHERE length(Sequence) >6  
>>       AND length(Sequence) <= 30
>>       GROUP BY Sequence 
>>       ORDER BY length(Sequence),Sequence
>
>    - **include query result column headers**: `Yes`
>
> 2. Click **Execute** and inspect the query results file after it turned green.<img src="../../../QT_output.png" width=60%> 

### Tabular to FASTA

> **Title column**: `1`
>
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

To continue processing this data, proceed to workflow 3 for Novel Peptide analysis.

## Part III


---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: Proteogenomics_NovelPeptide_Analysis
---

# Introduction

The third and the last proteogenomics workflow is for identifying the novel peptides using BlastP and to localize the peptides to its genomic coordinates. Inputs from both workflow 1 and 2 will be used in this workflow.

<img src="../../../images/Third_workflow.png" width=100%>

### Agenda
>
> In this tutorial, we will deal with:
>

> - _Inputs required_ 
>
> - _Interactive visualization of the Peptides_
>
> - _Classification of Novel Peptides_
>
> - _Summary of Identified Novel peptides_
>
>>      The inputs for this workflow are:
>>            - Tabular file – “Peptides from BlastP analysis”
>>            - Tabular file – “PeptideShaker_PSM”
>>            - Mz to sqlite
>>            - Genomic mapping sqlite

> All the files to run this workflow can be obtained from the second workflow output.Once the tabular output is created, 
> we convert this tabular report into a FASTA file. This can be achieved by using the Tabular to FASTA convertion tool.


Once Blast-P search is performed, it provides with a tabular output containing “Novel peptides”. Now this output is further processed by comparing the Novel Peptide output with the PSM report for selecting only distinct peptides which pass these steps.

# Analysis

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
> 2. Click **Execute** and inspect the query results file after it turned green. 
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
> 2. Click **Execute** and inspect the query results file after it turned green. 
>

### MVP

> The spectra belonging to these novel peptides can be viewed using MVP,this can be achieved by selecting the output from the mz to sqlite tool. 
> Here is a step by step proteogenomic view of the novel peptides obtained from running this workflow:

> 1) Click on the “Visualize in MVP application”, it will open up a new window for MVP.
> <img src="../../../images/mz2sqlite.png" width=20%>
>
> 2) Click on “Load from Galaxy”.
> <img src="../../../images/load_from_Galaxy.png" width=50%>
>
> 3) Select “Novel Peptides” from the right hand side.
> <img src="../../../images/novel_peptides_view.png" width=50%>
>
> 4) Select any peptide, For eg: ESSREALVEPTSESPRPALAR, and then click on “Selected Peptide PSMs”.
> <img src="../../../images/select_pep_PSM.png" width=50%>
>
> 5) If you scroll down, the PSM associated with the peptide will be displayed. By clicking on the PSM, the lorikeet 
> values will be shown. The lorikeet visualization is interactive, i.e the user can change the values or select any 
> parameter and click on Update button to view these changes.
>
> <img src="../../../images/Psm.png" width=40%>
>
> <img src="../../../images/lorikeet.png" width=70%>
>
> 6) For a Protein centric view, click on “View in Protein” , it will open up all the proteins associate with the 
> peptides. For eg: Select the “ESSREALVEPTSESPRPALAR” peptide and click on the first protein. The chromosome location 
> of the peptide will be displayed.
>
> <img src="../../../images/view_in_prot.png" width=30%>
>
> Once you click on protein it will show the list of proteins the belongs to the peptides.
>
> <img src="../../../images/select_protein.png" width=60%>
>
> Once you select the protein that you want to visualize you can click on the protein view.
>
> <img src="../../../images/PRoteinview.png" width=50%>
>
> 7) Clicking on the arrow marks will open up the IGV(js) visualization tool, where-in the genomic localization of the 
> peptide will be displayed.
>
> <img src="../../../images/select_IGV.png" width=80%>
>
> 8) To add tracks to your IGV viewer, click on “Add Track”. This will open up a list of tracks that are compatible 
> to view in your IGV viewer. For eg. Select the “Pep_gen_coordinate.bed” file and then click on “Load Track”.
> This will open up the bed will below the nucleotide sequence.
> <img src="../../../images/track_load.png" width=40%>
>
> 9) By clicking the wheel, you can select the “three frame translate” which will show the three frame translated 
> region of your sequence.
> <img src="../../../images/IGV_viewer.png" width=40%>
>
> 10) The IGV is inbuilt in the MVP viewer and is very interactive, you could also load more tracks such as the aligned 
> Bam file (from HISAT) or the identified pro bam file (one of the input file).
> MVP has many useful features beyond those covered in this workshop and is under active development.
> <img src="../../../images/tracks_align.png" width=70%>

The next tool in the workflow is the Peptide genomic coordinate tool which takes the novel peptides as the input along with the mztosqlite file and the genomic mapping sqlite file (obtained during creation of the database). This tool helps create a bed file with the genomic coordinate information of the peptides based on the sqlite files. 

### Peptide genomic Coordinate
Gets genomic coordinate of peptides based on the information in mzsqlite and genomic mapping sqlite files. This program 
loads two sqlite databases (mzsqlite and genomic mapping sqlite files) and calculates the genomic coordinates of the
peptides provided as input. This outputs bed file for peptides.
>
> 1. Peptide genomic Coordinate
>       - **Input**: `Peptide list file`, `mzsqlite sqlite DB file`, and `genomic mapping sqlite DB file` 
>       - **Output**: `Tabular BED file with all the columns`
> <img src="../../../images/pep_gen_cor.png" width=100%>
>
> mzsqlite file from: https://toolshed.g2.bx.psu.edu/repos/galaxyp/mz_to_sqlite/mz_to_sqlite/2.0.0 
> genome mapping sqlite file from: https://toolshed.g2.bx.psu.edu/view/galaxyp/translate_bed/038ecf54cbec
>
> <img src="../../../images/Output_PGC.png" width=50%>
>  2. Click **Execute** and inspect the query results file after it turned green. 
>

### Peppointer

Given chromosomal locations of peptides in a BED file, PepPointer classifies them as CDS, UTR, exon, intron, or intergene.

> 1. Peppointer
>      - **Choose the source of the GTF file** - `Locally Installed`
>              - **GTF file with the genome of interest** - `Mus_Musculus_GRCm38.90_Ensembl_GTF`
>      - **Input** - `Bed file from Peptide genomic coordinate tool`
> <img src="../../../images/Peppointer.png" width=80%>
>  This tool provides a bed output with the classification of the genomic location of the peptides.
> <img src="../../../images/Output_PP.png" width=50%>
> 2. Click **Execute** and inspect the query results file after it turned green. 

The final tool for this workflow is creating a tabular output that includes all the information that you get after running these workflows. The final summary output consists of the Peptide sequence, the spectra associated with the peptides, the protein accession number, chromosome number, Start and Stop of the genomic coordinate, the annotation, the genomic coordinate entry for viewing in Integrative Genomics Viewer (IGV), MVP or UCSC genome browser and the URL for viewing it on UCSC genome browser. This summary is created with the help of the query tabular tool.


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
<img src="../../../images/final_summary.png" width=80%>

This completes the proteogenomics workflow analysis. This training workflow uses mouse data but for any other organism the data and the workflow has to be modified accordingly.

This workflow is also available at z.umn.edu/proteogenomicsgateway.

This workflow was developed by the Galaxy-P team at the University of Minnesota.
For more information about Galaxy-P or our ongoing work, please visit us at www.galaxyp.org



