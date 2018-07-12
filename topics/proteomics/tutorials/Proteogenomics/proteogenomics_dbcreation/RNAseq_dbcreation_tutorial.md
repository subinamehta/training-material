---
layout: tutorial_hands_on
topic_name: proteomics
tutorial_name: Proteogenomics_RNAseq_db_creation
---

# Introduction

Proteogenomics is a combination of proteomics, genomics and transcriptomics data to identify peptides and to understand protein level evidence of gene expression. In this Proteogenomics tutorial we will create a Protein FASTA database using RNA sequencing files (FASTQ) and then perform Database searching of the created FASTA file with MS/MS data to identify Novel Peptides. We will then assign the genomic coordinate and annotation for these novel peptides as well as perform visualization of the data. 

Proteogenomics most commonly integrates RNA-Seq data, for generating customized protein sequence databases, with mass spectrometry-based proteomics data, which are matched to these databases to identify novel protein sequence variants. (Cancer Res. (2017); 77(21):e43-e46. doi: 10.1158/0008-5472.CAN-17-0331.)

In this tutorial, the proteins and the total RNA were obtained from the early development of B-cells from mouse. It was obtained at two developmental stages of B-cells, Ebf -/- pre-pro-B and Rag2 -/- pro-B. Please refer to the original study for details [Heydarian, M. et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4276347/).

### Agenda
>
> In this tutorial, we will deal with:
>
> 1. TOC
> {:toc}
!-- [[TOC]] --
> {: .agenda}

# Pretreatments

## Data upload

There are a many ways how you can upload your data. Three among these are:

*   Upload the files from your computer
*   Using a direct link
*   Import from the data library if your instance provides the files

In this tutorial, we will get the data from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.839701.svg)](https://doi.org/10.5281/zenodo.839701).

> ### {% icon hands_on %} Hands-on: Data upload and organization
>
> 1. Create a new history and name it something meaningful (e.g. *Proteogenomics tutorial*)
> 2. Import the three MGF MS/MS files and the FASTA sequence file from Zenodo.
>
>    > ### {% icon tip %} Tip: Importing data via links
>    >
>    > * Copy the link location
>    > * Open the Galaxy Upload Manager
>    > * Select **Paste/Fetch Data**
>    > * Paste the link into the text field. You can add multiple links, each on a separate line.
>    > * Press **Start**
>    {: .tip}
>
>    As default, Galaxy takes the link as name.
>
>    > ### {% icon comment %} Comments
>    > - Rename the datasets to a more descriptive name
>    {: .comment}
>
{: .hands_on}

# Analysis

The first workflow focuses on creating a FASTA Database created from RNA-seq data. There are two outputs from this workflow, a sequence database consisting of variants and known reference sequences and mapping files containing genomic and variant mapping data.

## Aligning FASTQ files on the human genome

The first tool in the workflow is the [HISAT2](http://ccb.jhu.edu/software/hisat) alignment tool. It maps next generation sequence reads to the reference genome. For running the HISAT2 tools there are two input files a RNA-seq file (.FASTQ) and a reference genome (GTF file format). The .gtf (Gene Transfer Format (GTF)) file is obtained from the Ensembl database.
This tool creates a .bam file.

#### HISAT2

> ### {% icon hands_on %} Hands-on: HISAT2
>
> 1. **HISAT2** {% icon tool %}: Run **HISAT2** with:
>    - **Source for the reference genome**: `Use a built-in genome` mm10
>    - **Single-end or paired-end reads**: `Single end` 
>    - **Input FASTQ files**: `FASTQ_ProB.fastqsanger`
>    - **Specify strand information**: `Unstranded`
>    > ### {% icon tip %} Tip: If you have paired inputs, select the paired end reads.
>   
>
>    Section **Summary Options**:
>
>    - Select `default parameters`

>   Section **Advanced Options**:

    - Select `default parameters`
    
>    > ### {% icon comment %} Comment
>    >
>    > Note that if your reads are from a stranded library, you need to choose the appropriate setting under Specify strand information above. For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF. With this option being used, every read alignment will have an XS attribute tag: '+' means a read belongs to a transcript on '+' strand of genome. '-' means a read belongs to a transcript on '-' strand of genome. (TopHat has a similar option, --library-type option, where fr - firststrand corresponds to R and RF; fr - secondstrand corresponds to F and FR.)

>    {: .comment}
>
> 2. Click **Execute**.
>
{: .hands_on}


#### FreeBayes

[FreeBayes]( https://github.com/ekg/freebayes) FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment.

> ### {% icon comment %} Comment
> Provided some BAM dataset(s) and a reference sequence, FreeBayes will produce a VCF dataset describing SNPs, indels, and complex variants in samples in the input alignments.
> By default, FreeBayes will consider variants supported by at least 2 observations in a single sample (-C) and also by at least 20% of the reads from a single sample (-F). These settings are suitable to low to high depth sequencing in haploid and diploid samples, but users working with polyploid or pooled samples may wish to adjust them depending on the characteristics of their sequencing data.
> FreeBayes is capable of calling variant haplotypes shorter than a read length where multiple polymorphisms segregate on the same read. The maximum distance between polymorphisms phased in this way is determined by the --max-complex-gap, which defaults to 3bp. In practice, this can comfortably be set to half the read length.
> Ploidy may be set to any level (-p), but by default all samples are assumed to be diploid. FreeBayes can model per-sample and per-region variation in copy-number (-A) using a copy-number variation map.
> FreeBayes can act as a frequency-based pooled caller and describe variants and haplotypes in terms of observation frequency rather than called genotypes. To do so, use --pooled-continuous and set input filters to a suitable level. Allele observation counts will be described by AO and RO fields in the VCF output.
{: .comment}

### {% icon hands_on %} Hands-on: Freebayes
>
> 1. **FreeBayes** {% icon tool %}:
>   - **Choose the source for the reference genome**: `Locally cached file`
>       - **Run in batch mode?**: `Run Individually`
>   - **BAM dataset**: `HISAT_Output.BAM`
>   - **Using reference genome**: `Human Dec.2013 (GRCh38/hg38)(hg38)`
>   - **Limit variant calling to a set of regions?**: `Do not Limit`
>   - **Choose parameter selection level**: `Simple diploid calling`
>   
>
### {% icon comment %} Comment
>  
>   Galaxy allows five levels of control over FreeBayes options, provided by the Choose parameter selection level menu option. These are:

> 1.Simple diploid calling: The simplest possible FreeBayes application. Equivalent to using FreeBayes with only a BAM input and no other parameter options.
> 2.Simple diploid calling with filtering and coverage: Same as #1 plus two additional options: -0 (standard filters: --min-mapping-quality 30 --min-base-quality 20 --min-supporting-allele-qsum 0 --genotype-variant-threshold 0) and --min-coverage.
> 3.Frequency-based pooled calling: This is equivalent to using FreeBayes with the following options: --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic. This is the best choice for calling variants in mixtures such as viral, bacterial, or organellar genomes.
> 4.Frequency-based pooled calling with filtering and coverage: Same as #3 but adds -0 and --min-coverage like in #2.
> Complete list of all options: Gives you full control by exposing all FreeBayes options as Galaxy parameters.
>      {: .comment}
>
> 2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../../../../view_icon.png)
>
{: .hands_on}

#### CustomProDB

[CustomProDB]( http://dx.doi.org/10.1093/bioinformatics/btt543) Generate custom protein FASTAs from exosome or transcriptome data.
The reference protein set can be filtered by transcript expression level (RPKM calculated from a BAM file), and variant protein forms can be predicted based on variant calls (SNPs and INDELs reported in a VCF file).

> ### {% icon comment %} Comment
> Annotations CustomProDB depends on a set of annotation files (in RData format) to create reference and variant protein sequences. Galaxy administrators can use the customProDB data manager to create these annotations to make them available for users.
> 
{: .comment}

### {% icon hands_on %} Hands-on: CustomProDB Generate protein FASTAs from exosome or transcriptome data
>
> 1. **CustomProDB** {% icon tool %}:
>   - **Will you select a genome annotation from your history or use a built-in annotation?**: `Use built in genome annotation`
>   - **Using reference genome**: `Human Ensembl 89 hsapiens (hg38/GRCh38.p10)`
>   - **BAM file**: `HISAT_Output.BAM`
>   - **VCF file**: `Freebayes.vcf`
>   - **Annotate SNPs with rsid from dbSNP**: `No`
>   - **Annotate somatic SNPs from COSMIC (human only)**: `No`
>   - **Transcript Expression Cutoff (RPKM)**: `1`
>   - **Create a variant FASTA for short insertions and deletions**: `Yes`
>   - **Create SQLite files for mapping proteins to genome and summarizing variant proteins**: `Yes`
>   - **Create RData file of variant protein coding sequences**: `Yes`


>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button]
>   
>
>       > ### {% icon comment %} Comment
>       Three FASTA files are created through the Custom ProDB tool, a variant FASTA file for short indels, a Single Amino acid Variant (SAV) FASTA file, an Sqlite file (genome mapping and variant mapping) for mapping proteins to genome and a RData file for variant protein coding sequences.

Similar to the genomic mapping, a variant mapping file is also created from CustomProDB. This sqlite file is also converted to tabular and made SearchGUI compatible. This variant annotation file will be used to visualize the variants in the Multi-omics visualization Platform (in-house visualization platform developed by Galaxy-P senior developers).


#### StringTie

[StringTie](http://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. 
> Its input can include not only the alignments of raw reads used by other transcript assemblers, but also alignments of longer sequences that have been assembled from those reads. In order to identify differentially expressed genes between experiments, StringTie's output can be processed by specialized software like Ballgown, Cuffdiff or other programs (DESeq2, edgeR, etc.).


1. **StringTie** {% icon tool %}:
> **StringTie transcript assembly and quantification**
>   - **Input mapped reads**: `HISAT_Output.BAM`
>   - **Specify strand information**: `Unstranded`
>   - **Use a reference file to guide assembly?**: `Use Reference GTF/GFF3`
>   - **Reference file**: `Use file from History`
>       - **GTF/GFF3 dataset to guide assembly**: `Homo_sapiens.GRCh38.83.gtf`
>   - **Use Reference transcripts only?**: `No`
>   - **Output files for differential expression?**: `No additional output`
>   - **Output coverage file?**: `No`
>   - **Advanced Options**: `Default Parameters`


>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../images/view_data_icon.png)

> ### {% icon comment %} Comment
>       StringTie takes as input a BAM (or SAM) file of paired-end RNA-seq reads, which must be sorted by genomic location (coordinate position). This file contains spliced read alignments and can be produced directly by programs such as HISAT2. We recommend using HISAT2 as it is a fast and accurate alignment program. Every spliced read alignment (i.e. an alignment across at least one junction) in the input BAM file must contain the tag XS to indicate the genomic strand that produced the RNA from which the read was sequenced. Alignments produced by HISAT2 (when run with --dta option) already include this tag, but if you use a different read mapper you should check that this XS tag is included for spliced alignments.

> NOTE: be sure to run HISAT2 with the --dta option for alignment (under Spliced alignment options), or your results will suffer.

> Also note that if your reads are from a stranded library, you need to choose the appropriate setting under Specify strand information above. As, if Forward (FR) is selected, StringTie will assume the reads are from a --fr library, while if Reverse (RF) is selected, StringTie will assume the reads are from a --rf library, otherwise it is assumed that the reads are from an unstranded library (The widely-used, although now deprecated, TopHat had a similar --library-type option, where fr-firststrand corresponded to RF; fr-secondstrand corresponded to FR). If you don't know whether your reads are from are a stranded library or not, you could use the tool RSeQC Infer Experiment to try to determine.

> As an option, a reference annotation file in GTF/GFF3 format can be provided to StringTie. In this case, StringTie will prefer to use these "known" genes from the annotation file, and for the ones that are expressed it will compute coverage, TPM and FPKM values. It will also produce additional transcripts to account for RNA-seq data that aren't covered by (or explained by) the annotation. Note that if option -e is not used the reference transcripts need to be fully covered by reads in order to be included in StringTie's output. In that case, other transcripts assembled from the data by StringTie and not present in the reference file will be printed as well.

NOTE: we highly recommend that you provide annotation if you are analyzing a genome that is well-annotated, such as human, mouse, or other model organisms.

#### GffCompare
> [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) compare and evaluate the accuracy of RNA-Seq transcript assemblers (Cufflinks, Stringtie). * collapse (merge) duplicate transcripts from multiple GTF/GFF3 files (e.g. resulted from assembly of different samples) * classify transcripts from one or multiple GTF/GFF3 files as they relate to reference transcripts provided in a annotation file (also in GTF/GFF3 format)

> The original form of this program is also distributed as part of the Cufflinks suite, under the name ["CuffCompare"] (http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/). Most of the options and parameters of CuffCompare are supported by GffCompare, while new features will likely be added to GffCompare in the future.

1. **GffCompare compare assembled transcripts to a reference annotation** {% icon tool %}:
>   -**GTF inputs for comparison**`Stringtie_outut.gtf`
>   - **Use Reference Annotation**: `Homo_sapiens.GRCh38.83.gtf`
>   - **Reference Annotation**: `Unstranded`

>   - **Ignore reference transcripts that are not overlapped by any input transfrags**: `No`
>   - **Ignore input transcripts that are not overlapped by any reference transcripts**: `No`
>   - **Use Sequence Data**: `No`
>   - **discard (ignore) single-exon transcripts**: `No`
>   - **Max. Distance for assessing exon accuracy**: `100`
>   - **Max distance for transcript grouping**: `100`
>   - **discard intron-redundant transfrags sharing 5'**: `No`

>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../images/view_data_icon.png)

>  ### {% icon comment %} Comment
>       A notable difference from GffCompare is that when a single query GTF/GFF file is given as input, along with a reference annotation (-r option), gffcompare switches into "annotation mode" and it generates a .annotated.gtf file instead of the .merged.gtf produced by CuffCompare with the same parameters. This file has the same general format as CuffCompare's .merged.gtf file (with "class codes" assigned to transcripts as per their relationship with the matching/overlapping reference transcript), but the original transcript IDs are preserved, so gffcompare can thus be used as a simple way of annotating a set of transcripts.

> Another important difference is that the input transcripts are no longer discarded when they are found to be "intron redundant", i.e. contained within other, longer isoforms. CuffCompare had the -G option to prevent collapsing of such intron redundant isoforms into their longer "containers", but GffCompare has made this the default mode of operation (hence the -G option is no longer needed and is simply ignored when given).

#### Convert gffCompare annotated GTF to BED for StringTie results

> Convert a GFFCompare annotated GTF file to BED format.

1. **GffCompare compare assembled transcripts to a reference annotation** {% icon tool %}:
> **GTF annotated by gffCompare**
> **filter gffCompare class_codes to convert** 
> `j : Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript
> e : Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.
> i : A transfrag falling entirely within a reference intron
> p : Possible polymerase run-on fragment (within 2Kbases of a reference transcript)
> u : Unknown, intergenic transcript`

>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../images/view_data_icon.png)

#### Translate BED transcripts cDNA in 3frames or CDS
> Translate transcripts from the input BED file into protein sequences.

> The genomic sequence:
1.may be supplied in an extra column in the BED input file
2.retrieved from a twobit genomic reference file
3.retrieved from the Ensembl REST API for Ensembl transcripts

1. **Translate BED transcripts cDNA in 3frames or CDS** {% icon tool %}:
**A BED file with 12 columns**: `Convert gffCompare annotated GTF to BED`
**Source for Genomic Sequence Data** `Locally cached File`
**Select reference 2bit file** `hg38`
**BED Filtering Options** `default`
**Translation Options ** `default`
**FASTA ID Options** `default`

>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../images/view_data_icon.png)

> ### {% icon comment %} Comment
> **INPUTS**

> BED file with at least the standard 12 columns
> Genome reference in twobit format (optional)

> **OUTPUTS**

>FASTA of transcript translations
>BED with the genomic location of the translated protein. The added 13th column contains the protein sequence.

> **OPTIONS**

> Feature translation
    - cDNA - three frame translations of the cDNA sequences with an output for each sequence between STOP codons
    - CDS - three frame translations of CDS (coding sequence defined by thickStart and thickEnd in the BED file)
> Translation filtering
    - can be trimmed to a Methionine start codon
    - can be split into peptides by an enzyme digestion
    - must exceed specified minimum length
> BED Filtering
    - genomic regions
    - ensembl biotype if the BED contains the 20 columns as retrieved from the Ensembl REST API
    
#### bed to protein map genomic location of proteins for MVP

Convert a BED format file of the proteins from a proteomics search database into a tabular format for the Multiomics Visualization Platform (MVP).

>   1. **A BED file with 12 columns, thickStart and thickEnd define protein coding region**: `Translate cDNA_minus_CDS`

>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../images/view_data_icon.png)

> ### {% icon comment %} Comment

> The tabular output can be converted to a sqlite database using the Query_Tabular tool.

> The sqlite table should be named: feature_cds_map The names for the columns should be: name,chrom,start,end,strand,cds_start,cds_end

> This SQL query will return the genomic location for a peptide sequence in a protein (multiply the animo acid position by 3 for the cds location):


### Creating FASTA Database:

The Protein database downloader tool is used to download the FASTA database from UNIPROT and cRAP database containing known/reference mouse proteins.

#### FASTA Merge Files and Filter Unique Sequences Concatenate FASTA database files together

> Concatenate FASTA database files together.
> If the uniqueness criterion is "Accession and Sequence", only the first appearence of each unique sequence will appear in the output. Otherwise, duplicate sequences are allowed, but only the first appearance of each accession will appear in the output.
> The default accession parser will treat everything in the header before the first space as the accession.

> 1. **Run in batch mode?**: `Merge individual FASTAs (output collection if input is collection)`

> *Input FASTA File(s)* : ` Input Custom ProDB Fasta File output`
                         > 1.HISAT_Output.rpkm
                         > 2.HISAT_Output.snv
                         > 3.HISAT_Output.indel

**How are sequences judged to be unique?**:`Accession and Sequence`
**Accession Parsing Regex**: `^>([^ |]+).*$`

>   2. Click **Execute** and inspect the resulting files after they turned green with the **View data** icon:
>     ![View data button](../../../images/view_data_icon.png)

> ### {% icon comment %} Comment
The regex text manipulation tool is used to manipulate the FASTA file to make it searchGUI compatible. The “FASTA Merge Files and Filter Unique Sequences Concatenate FASTA databases together” tool is used to merge the databases obtained from the CustomProDB and translate Bed tool along with the Uniprot and cRAP databases.


#### Recieving the list of peptides: Query Tabular

In order to use *Unipept*, a list containing the peptide sequences has to be generated.
The tool **Query Tabular** can load tabular data (the PSM report in this case) into a SQLite data base.
As a tabular file is being read, line filters may be applied and an SQL query can be performed.

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`:
>    - **Tabular Dataset for Table**: The PSM report
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `^\d`
>        - **action for regex match**: `include line on pattern match`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`
>    - **Specify Column Names (comma-separated list)**: `id,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>
>        > ### {% icon comment %} Comment
>        >
>        > By default, table columns will be named: c1,c2,c3,...,cn (column names for a table must be unique).
>        > You can override the default names by entering a comma separated list of names, e.g. `,name1,,,name2`
>        > would rename the second and fifth columns.
>        >
>        > Check your input file to find the settings which best fits your needs.
>        {: .comment}
>
>    - **Only load the columns you have named into database**: `Yes`
>
>    - **Save the sqlite database in your history**: `Yes`
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
>          SELECT distinct sequence
>
>          FROM psm
>
>          WHERE confidence >= 95
>
>          ORDER BY sequence
>
>    > ### {% icon question %} Questions
>    >
>    > The SQL query might look confusing at first, but having a closer look should clarify a lot.
>    >
>    > 1. What does `FROM psm` mean?
>    > 2. What need to be changed if we only want peptides with a confidence higher then 98%?
>    >
>    >    > ### {% icon solution %} Solution
>    >    > 1. We want to read from table "psm". We defined the name before in the "Specify Name for Table" option.
>    >    > 2. We need to change the value in line 3: "WHERE validation IS NOT 'Confident' AND confidence >= 98"
>    >    {: .solution }
>    {: .question}
>
>    - **include query result column headers**: `No`
>
> 2. Click **Execute** and inspect the query results file after it turned green. If everything went well, it should look similiar:
>
>     ![Query Tabular output showing the peptides](../../../images/query_tabular_1.png "Query Tabular output")
>
{: .hands_on}

While we can proceed with this list of peptides, let's practice using the created SQLite database for further queries.
We might not only be interested in all the distinct peptides, but also on how many PSMs a single peptide had.
Therefore we can search the database for the peptides and count the occurrence without configuring the tables and columns again:

> ### {% icon hands_on %} Hands-on: SQLite to tabular
>
> 1. **SQLite to tabular** {% icon tool %}: Run **SQLite to tabular** with:
>
>    - **SQL Query**:
>
>          SELECT sequence as "peptide", count(id) as "PSMs"
>
>          FROM psm
>
>          WHERE confidence >= 95
>
>          GROUP BY sequence
>
>          ORDER BY sequence
>
> 2. Click **Execute**. The resulting file should have two columns, one with the distinct peptides, the other with the count number of PSMs.
>
{: .hands_on}


#### Retrieve taxonomy for peptides: Unipept

The generated list of peptides can now be used to search via *Unipept*.
We do a taxonomy analysis using the UniPept pept2lca function to return the taxonomic lowest common ancestor for each peptide:

> ### {% icon hands_on %} Hands-on: Unipept
>
> 1. **Unipept** {% icon tool %}: Run **Unipept** with:
>
>    - **Unipept application**: `pept2lca: lowest common ancestor`
>    - **Peptides input format**: `tabular`
>    - **Tabular Input Containing Peptide column**: The query results file.
>    - **Select column with peptides**: `Column 1`
>    - **Choose outputs**: Select `tabular` and `JSON taxonomy tree`
>
> 2. Click **Execute**. The history should grow by two files. View each to see the difference.
>
>       > ### {% icon comment %} Comment
>       >
>       > The JSON (JavaScript Object Notation) file contains the same information as the tabular file but is not comfortably human readable.
>       > Instead, we can use it to use JavaScript libraries to visualize this data.
>       {: .comment}
>
> 3. Visualize the data:
>
>    - Click on the JSON output file from the *Unipept* tool to expand it. Click on the **Visualize** button and select **Unipept Tree viewer**:
>
>       ![Visualize button](../../../images/visualize_button.png)
>
>    - A new window should appear with a visualization of the taxonomy tree of your data. Use the mouse wheel to scroll in and out and click on nodes to expand or collapse them:
>
>       ![Unipept Tree viewer visual output](../../../images/unipept_tree_viewer.png "Interactive visualization from the Unipept Tree viever plugin")
>
{: .hands_on}

## Genus taxonomy level summary

The tabular *Unipept* output lists the taxonomy assignments for each peptide. To create a meaningful summary, the **Query Tabular** tool is
once again used, aggregating the number of peptides and PSMs for each genus level taxonomy assignment:

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The PSM report
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `^\d`
>        - **action for regex match**: `include line on pattern match`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `psm`
>    - **Specify Column Names (comma-separated list)**: `,,sequence,,,,,,,,,,,,,,,,,,,,confidence,validation`
>
>    - **Only load the columns you have named into database**: `Yes`
>
> 2. Repeat this step to have a second **Database Table**:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The **Unipept** `tabular`/`tsv` output
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `#peptide`
>        - **action for regex match**: `exclude line on pattern match`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `lca`
>    - **Specify Column Names (comma-separated list)**: `peptide,,,,,,,,,,,,,,,,,,,,,genus`
>
>    - **Only load the columns you have named into database**: `Yes`
>
>    - **Save the sqlite database in your history**: `Yes`
>
>    - **SQL Query to generate tabular output**:
>
>          SELECT lca.genus,count(psm.sequence) as "PSMs",count(distinct psm.sequence) as "DISTINCT PEPTIDES"
>
>          FROM psm LEFT JOIN lca ON psm.sequence = lca.peptide
>
>          WHERE confidence >= 95
>
>          GROUP BY lca.genus
>
>          ORDER BY PSMs desc, 'DISTINCT PEPTIDES' desc
>
>
> 2. Click **Execute** and inspect the query results file after it turned green:
>
>     ![Query Tabular output showing gene, PSMs and distinct peptides](../../../images/metaproteomics_summary.png "Query Tabular output")
>
{: .hands_on}

## Functional Analysis

Recent advances in microbiome research indicate that functional characterization via metaproteomics analysis has the potential to accurately
measure the microbial response to perturbations. In particular, metaproteomics enables the estimation of the function of the microbial
community based on expressed microbial proteome.

In the following chapter, a functional analysis will be performed using the **UniPept** application `pept2prot` in order to match the list of peptides with the correlated Gene Ontology terms.
This allows to get an insight of the **biological process**, the **molecular function** and the **cellular component** related to the sample data.

> ### {% icon comment %} Gene Ontology Consortium
>
> The [Gene Ontology Consortium](http://www.geneontology.org/) provides with its Ontology a framework for the model of biology.
> The GO defines concepts/classes used to describe gene function, and relationships between these concepts. It classifies functions along three aspects:
>
>
> - **molecular function**
>
>   - molecular activities of gene products
>
> - **cellular component**
>
>   - where gene products are active
>
> - **biological process**
>
>   - pathways and larger processes made up of the activities of multiple gene products.
>
> [more information](http://geneontology.org/page/ontology-documentation)
>
{: .comment}

#### Data upload

For this tutorial, a tabular file containing the relevant GO terms has been created. It contains the GO aspect, the ID and the name.
It is available at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.839701.svg)](https://doi.org/10.5281/zenodo.839701).

> ### {% icon hands_on %} Hands-on: Data upload
>
> 1. Import the file `Gene_Ontology_Terms.tabular` from Zenodo.
>
>    > ### {% icon tip %} Tip: Setting file metadata on upload
>    >
>    > In the upload window of Galaxy you can set the filetype and related genome of the file you're uploading in the corresponding columns beforehand.
>    > This might be handy if the automatic detection of the filetype didn't work out perfectly or if you want to avoid setting the genome later on, especially for multiple files.
>    >
>    {: .tip}
>
>    As default, Galaxy takes the link as name.
>
>    > ### {% icon comment %} Comments
>    > - Rename the datasets to a more descriptive name, e.g. `Gene Ontology Terms`
>    {: .comment}
>
>
{: .hands_on}

> ### {% icon tip %} Tip: Creating your own Gene Ontology list
>
> The latest Gene Ontology can be downloaded [here](http://geneontology.org/page/download-ontology) as a text file in the `OBO` format.
> `OBO` files are human-readable (in addition to machine-readable) and can be opened in any text editor. They contain more information than just the name and aspect.
>
> In order to receive a file like we use in the tutorial for your own analysis, different tools are available to extract information from `OBO` files,
> one of them being [ONTO-PERL](https://doi.org/10.1093/bioinformatics/btn042).
> An example file with all GO terms from 08.07.2017 named `Gene_Ontology_Terms_full_07.08.2017.tabular` can be found on the [Zenodo repository](https://doi.org/10.5281/zenodo.839701) of this tutorial as well.
>
{: .tip}

#### Retrieve GO IDs for peptides: Unipept

The **UniPept** application `pept2prot` can be used to return the list of proteins containing each peptide.
The option `retrieve extra information` option is set to `yes` so that we retrieve Gene Ontology assignments (`go_references`)
for each protein.

> ### {% icon hands_on %} Hands-on: Unipept
>
> 1. **Unipept** {% icon tool %}: Run **Unipept** with:
>
>    - **Unipept application**: `pept2prot: UniProt entries containing a given tryptic peptide`
>    - **retrieve extra information**: `Yes`
>    - **Peptides input format**: `tabular`
>    - **Tabular Input Containing Peptide column**: The first query results file.
>    - **Select column with peptides**: `Column 1`
>    - **Choose outputs**: Select `tabular`
>
> 2. Click **Execute**.
>
> 3. inspect the result:
>
>    - The output should be a tabular file containing a column labeled `go_references`. This is what we're looking for.
>
{: .hands_on}


#### Combine all information to quantify the GO results

As a final step we will use **Query Tabular** in a more sophisticated way to combine all information to quantify the GO analysis. The three used file and the extracted information are:

- **Gene Ontology Terms**:
    - `go_id` to match with **Normalized UniPept output**
    - The GO `aspect` to group the results in three separate files
    - The GO `description` to annotate the results
- **Normalized UniPept output**:
    - `peptide` to match with **PSM Report** and to count distinct peptides per GO term
    - `go_reference` to match with **Gene Ontology Terms**
- **PSM Report**:
    - `sequence` to match with **Normalized UniPept output**
    - `id` to count distinct PSM's per GO term

> ### {% icon hands_on %} Hands-on: Query Tabular
>
> 1. **Query Tabular** {% icon tool %}: Run **Query Tabular** with:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The `Gene Ontology Terms` file
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `skip leading lines`
>        - **Skip lines**: `1`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `go`
>    - **Specify Column Names (comma-separated list)**: `aspect,go_id,description`
>    - **Table Index**: Click on `+ Insert Table Index`:
>        - **This is a unique index**: `No`
>        - **Index on Columns**: `aspect,go_id`
>
>
> 2. Repeat this step to have a second **Database Table**:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The **Unipept** `tabluar`/`tsv` output
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `skip leading lines`
>        - **Skip lines**: `1`
>    - Add another Filter: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `prepend a line number column`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `bering_prot`
>    - **Specify Column Names (comma-separated list)**: `id,peptide,uniprot_id,taxon_id,taxon_name,ec_references,go_references,refseq_ids,refseq_protein_ids,insdc_ids,insdc_protein_ids`
>    - **Table Index**: Click on `+ Insert Table Index`:
>        - **This is a unique index**: `No`
>        - **Index on Columns**: `id,peptide`
>
> 3. Repeat this step to have another **Database Table**:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The same **Unipept** `tabluar`/`tsv` output
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `skip leading lines`
>        - **Skip lines**: leave blank
>    - Add another Filter: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `prepend a line number column`
>    - Add another Filter: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `select columns`
>        - **enter column numbers to keep**: `1,7`
>    - Add another Filter: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `normalize list columns, replicates row for each item in list`
>        - **enter column numbers to normalize**: `2`
>        - **List item delimiter in column**: ` ` (a single blank character)
>
>    > ### {% icon comment %} Comments
>    > - The UniPept result file can contain multiple GO IDs in a single row. In order to create a normalized table of this data, these rows will be split so each record contains only one GO ID.
>    {: .comment}
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `bering_prot_go`
>    - **Specify Column Names (comma-separated list)**: `id,go_reference`
>    - **Table Index**: Click on `+ Insert Table Index`:
>        - **This is a unique index**: `No`
>        - **Index on Columns**: `go_reference,id`
>
> 4. Repeat this step to have another **Database Table**:
>
>    - **Database Table**: Click on `+ Insert Database Table`
>    - **Tabular Dataset for Table**: The `PSM Report`
>
>    Section **Filter Dataset Input**:
>
>    - **Filter Tabular Input Lines**: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `by regex expression matching`
>        - **regex pattern**: `^\d`
>        - **action for regex match**: `include line on pattern match`
>    - Add another Filter: Click on `+ Insert Filter Tabular Input Lines`:
>    - **Filter By**: Select `select columns`
>        - **enter column numbers to keep**: `1,3,23,24`
>
>    Section **Table Options**:
>
>    - **Specify Name for Table**: `bering_psms`
>    - **Specify Column Names (comma-separated list)**: `id,sequence,confidence,validation`
>    - **Only load the columns you have named into database**: `Yes`
>    - **Table Index**: Click on `+ Insert Table Index`:
>        - **This is a unique index**: `No`
>        - **Index on Columns**: `sequence,id`
>
>    - **Save the sqlite database in your history**: `Yes`
>
>    - **SQL Query to generate tabular output**:
>
>          SELECT sequence as "peptide", count(id) as "PSMs"
>
>          FROM bering_psms
>
>          WHERE confidence >= 95
>
>          GROUP BY sequence
>
>          ORDER BY sequence
>
>
> 5. Click **Execute**.
>
{: .hands_on}

With this we have combined all the data into a single database which we can now query to extract the desired information with **SQLite to tabular**:

> ### {% icon hands_on %} Hands-on: SQLite to tabular
>
> 1. **SQLite to tabular** {% icon tool %}: Run **SQLite to tabular** with:
>
>    - **SQLite Database**: The created SQLite database from the former step
>    - **SQL Query**:
>
>          SELECT go.description, 
>
>          count(distinct bering_psms.sequence) as "bering_peptides", count(distinct bering_psms.id) as "bering_psms" 
>
>          FROM go JOIN bering_prot_go ON go.go_id = bering_prot_go.go_reference JOIN bering_prot on bering_prot_go.id = bering_prot.id JOIN 
>
>          bering_psms ON bering_prot.peptide = bering_psms.sequence
>
>          WHERE go.aspect = 'molecular_function'
>
>          GROUP BY go.description
>
>          ORDER BY  bering_peptides desc,bering_psms desc
>
> 2. Click **Execute**.
> 3. Repeat these steps two times by replacing `molecular_function` in the fifth row of the SQL query by `biological_process` and `cellular_component`.
>
{: .hands_on}

With these three resulting files the functional analysis of this tutorial is finished. Each record contains the name of a GO term, the amount of peptides related to it and the amount of PSMs for these peptides. 

> ### {% icon comment %} References
>
> - [Dataset](https://www.ncbi.nlm.nih.gov/pubmed/27824341) and [SixGill software](https://www.ncbi.nlm.nih.gov/pubmed/27396978)
>
> - [Galaxy workflows for metaproteomics](https://www.ncbi.nlm.nih.gov/pubmed/26058579)
>
> - [Metaproteomics community effort](https://z.umn.edu/gcc2017mporal)
>
> - [Unipept](https://www.ncbi.nlm.nih.gov/pubmed/28552653)
>
{: .comment}

