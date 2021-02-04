
# Comparative_Coverage_Analysis
Repository of custom python scripts used in "A transcriptomic pipeline adapted for genomic sequence discovery of germline restricted sequence in zebra finch, Taeniopygia guttata", Asalone et al. 2021.

## Installation

* Clone the repository.

* Dependencies:
	* Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* Samtools http://samtools.sourceforge.net/
	* Stringtie https://ccb.jhu.edu/software/stringtie/
	* R 3.0.1 https://cran.rstudio.com/
	* RStudio https://rstudio.com/products/rstudio/download/
	* R packages
	  * Ballgown
	  * ggplot2
	  * tidyverse
	* Python v. 2.7.16 https://www.python.org/downloads/release/python-2716/
	* Biopython https://biopython.org/
	
## Protocol 1. Read Mapping

### 1.1 Bowtie2 
To index the assembly file:
```{r, eval=FALSE}
$ bowtie2-build <assembly.fasta> <index_prefix>
```

Map raw reads to the assembly:
```{r, eval=FALSE}
$ bowtie2 -x <index_prefix> -U <raw_reads.fastq> -S <output.sam>
```

### 1.2 Samtools
Convert the sam to a bam file:
```{r, eval=FALSE}
$ samtools view -bu <output.sam> ><output.bam>
```

Coordinate sort the bam file:
```{r, eval=FALSE}
$ samtools sort -T temp -O BAM -o <coor_sort_output.bam> <output.bam>
```

## Protocol 2. Stringtie-Ballgown

### 2.1 Stringtie
Make the annotation file with python script
```{r, eval=FALSE}
$ /PATH_to_directory/buildgff_assembly.py <assembly.fasta>
```
Will output a file named <assembly.fasta>.gff  in which each contig is denoted as a transcript. 

Run Stringtie
```{r, eval=FALSE}
$ /PATH_to_Stringtie/stringtie -te -b /PATH/output_folder/ <coor_sort_output.bam> -G <assembly.fasta.gff> -o <output.gtf>
```
Note: be sure to name the output folder with a prefix then an identifier: i.e. SampleTestis1 or SampleLiver1

### 2.2 Adjust files
Stringtie will output the output.gtf along with five files within the output_folder.
The five files will be named:

1. e2t.ctab

2. e_data.ctab

3. i2t.ctab

4. i_data.ctab

5. t_data

Two of these files will need to be adjusted to tell Ballgown that there are no introns. These files can be edited on the command line using `nano` or simply opening the files with a text editor such as TextEdit.

i2t.ctab will originally look like:
```{r, eval=FALSE}
i_id    t_id
```
This file should be adjusted to look like:
```{r, eval=FALSE}
i_id    t_id 
0       0
```

i_data.ctab will originally look like:
```{r, eval=FALSE}
i_id    chr     strand  start   end     rcount  ucount  mrcount
```
This file should be adjusted to look like: 
```{r, eval=FALSE}
i_id    chr     strand  start   end     rcount  ucount  mrcount
0       0       .       0       0       0       0       0
```

### 2.3 Create Phenotype Data
Ballgown requires a phenotype data file that contains two columns:
```{r, eval=FALSE}
id    description
```
id = the output_folder names from all stringtie runs.
description = categories of samples (i.e. somatic and germline)

The file may look something like this:
```{r, eval=FALSE}
id               description
sampleTestis1    germline
sampleLiver1     somatic
```
This file can be created in excel, however, be sure to save it as .csv

### 2.4 Upload Phenotype Data into R studio
```{r, eval=FALSE}
pheno_data <- read.csv(file = "/PATH_to_file/pheno_data.csv")
```

### 2.5 Ballgown
```{r, eval=FALSE}
bg <- ballgown(dataDir = "/PATH_to_output_folders/",
               samplePattern = 'sample', 
               pData = pheno_data,
               meas = c('FPKM'))
```
If you would like to output the individual FPKMs for each contig for each sample:
```{r, eval=FALSE}
contig_FPKM = texpr(bg, 'FPKM')
write.csv(contig_FPKM, file = "/PATH_where_to_write/contig_FPKM.csv")
```

### 2.6 Stattest
To calculate the differential expression (or fold change):
```{r, eval=FALSE}
diffexp <- stattest(bg, feature = 'transcript', meas = 'FPKM', covariate = 'description', getFC = TRUE)
```

Check the fold change values to ensure that the correct denominator was used (i.e. germline FPKM/somatic FPKM). If the values seem to be backwards (i.e. somatic FPKM/germline FPKM) simply create a new column for inverted fold change:
```{r, eval=FALSE}
diffexp$invert_fc <- 1 / diffexp$fc
```

### 2.7 Adding Contig ID and Lengths
The data table will not contain the contig ID or lengths. These can be added by completing the following steps:
Read in one of the e_data.ctab files since it contains both the id used in Ballgown and contig id:
```{r, eval=FALSE}
contig_id <- read_tsv(file = "/PATH_to_output_folder/e_data.ctab")
```
Adjust contig_id data table:
```{r, eval=FALSE}
contig_id %>%
  mutate("id" = as.factor(e_id)) %>%
  select(c("id", "chr")) ->
  contig_id
```
Join contig_id with diffexp:
```{r, eval=FALSE}
left_join(diffexp, contig_id, by = "id") ->
  diffexp_edit
```

Get lengths of contigs:
```{r, eval=FALSE}
$ /PATH_to_script/contig_lengths.py <assembly.fasta>
```
Outputs a file titles <assembly.fasta>_lengths.csv

Read into RStudio:
```{r,eval=FALSE}
contig_length <- read.csv(file = "/PATH_to_file/assembly.fasta_lengths.csv")
```
Join contig_length with diffexp_edit:
```{r, eval=FALSE}
left_join(diffexp_edit, contig_length, by = "chr") ->
  diffexp_final
```


### 2.8 Graphing Volcano Plot
```{r, eval=FALSE}
ggplot(diffexp_final, 
       mapping = aes(x = log2(invert_fc), y = log10(qval)) +
  geom_point() +
  scale_y_reverse() +
  geom_vline(xintercept = log2(2)) +
  geom_hline(yintercept = log10(0.05))+
  theme_bw() + 
  ylab("log10(q-value)") +
  xlab("log2(fold change)")->
  vol_plot
vol_plot
ggsave("vol_FPKM", plot = vol_plot, device = "pdf", path = "/PATH_to_save/")
```
If you would like to change color, shapes, or alpha values for individual types of data (i.e. qPCR validated) use the following within the ggplot:
```{r, eval=FALSE}
  scale_color_manual(values = c(#fill in based on number of types)) +
  scale_shape_manual(values = c(#fill in based on number of types))+
  scale_alpha_manual(values = c(#fill in based on number of types)) +
  scale_size_manual(values = c(#fill in based on number of types))+
```

### 2.9 Filter for High Confidence Set
To filter for the high confidence set of outliers (i.e. GRC contigs):
```{r, eval=FALSE}
diffexp_final %>%
  filter(length > 2000) %>%
  filter(qval < 0.05) %>%
  filter(invert_fc > 2)->
  high_conf_outliers
```
May need to adjust the cut-offs based on datasets. 

Do download the high confidence set: 
```{r, eval=FALSE}
write.csv(high_conf_outliers, file = "high_conficence_outliers.csv")
```

## Supplemental Scripts

### buildgff.py
```{r, eval=FALSE}
$ /PATH_to_script/buildgff.py <blast_outfmt6.txt>
```
Builds a gff from a blast output in outfmt 6. 

## Scripts used for Reciprocal BLAST Match
### RBM.py
```{r, eval=FALSE}
$ /PATH_to_script/RBM.py <blast_output1> <blast_output2>
```
Outputs a file with a list of the reciprocal BLAST matches. This script does not require -max_target_seqs 1, instead it just takes the top hit from each blast output for each query. For reciprocal BLAST matches, both the query and subject will be printed. In cases where there is only a one-way match, the query and "one-way" is printed. If a query is not in a file a "0" is printed. 
### get_names_RBM.py
```{r, eval=FALSE}
$ /PATH_to_script/get_names_RBM.py <contig_list_to_get> <RBM.txt>
```
To pull out a set of RBM for contigs of interest, create a list of the desired contigs' RBM. The script is created to identify the contig from the query column, however, if the known contig IDs are in the subject column simply change line 17 from `if contigDict.has_key(q1):` to `if contigDict.has_key(s1):`
