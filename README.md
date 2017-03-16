BICA (bone-in-culture array)
============

This repository contains custom codes used for BICA paper, which has been accepted by **Nature Communications**.

* **SupplementaryCode1.sh**: bash code for reference index building, and read alignment and quantification.
* **boneRefSamples folder**
	* **Ref folder**: RNA-seq profiling of 13 kinds of bone marrow stroma cells from [E-MTAB-2923](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2923/);
	* **BICA_DESeq.txt**: gene signature file of 13 kinds of bone marrow stroma cells for [CIBERSORT](https://cibersort.stanford.edu/index.php);
	* **SupplementaryCode2.Rmd**: R markdown code for generating gene signature file using RNA-seq profiling files in **Ref folder**.
* **ssGSEAproj folder**
	* **SupplementaryCode3.Rmd**: R markdown code for generating the input file for [ssGSEA](https://www.ncbi.nlm.nih.gov/pubmed/19847166);
	* **ssGSEA output files**: output files using Gene Ontology Biological Process, Hallmarks, and Oncogenic gene sets database.
* XenomeCorrect.py: custom python script that correct the output **fastq** file from **Xenome**.

## Data Download

The RNA-seq data have be deposited to NCBI GEO database. The GEO accession number is [GSE84114](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84114).

There are two processed files. The **sra** compressed files are stored in short read archive (sra) database.

* [GSE84114\_BICA_RNAseq.RData.gz](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84114/suppl/GSE84114_BICA_RNAseq.RData.gz) is the read count table using [RSEM v1.2.28](http://deweylab.github.io/RSEM/). Since **RSEM** uses Expectation-Maximization algorithm, desimals may appear in the gene/isoform expression count table.

* [GSE84114\_BICA_RNAseq.Normalized.RData.gz](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84114/suppl/GSE84114_BICA_RNAseq.Normalized.RData.gz) is the normalized data of human cancer cells and mouse stroma cells. Two normalization approaches in [DESeq2 v1.14.1](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) were used: regularized log transformation (**rld**) or variance stabilizing transformation (**vsd**).

* Raw RNA-seq reads are saved in **.sra** format. You need install [sra-tools v2.3.5-2](https://github.com/ncbi/sra-tools) and run the command below to extract the pair-end sequencing data.

```
fastq-dump -I --split-3 SRRXXXXXXX.sra
```

## In silico sorting

We used [Xenome v1.0.1](https://academic.oup.com/bioinformatics/article/28/12/i172/269972/Xenome-a-tool-for-classifying-reads-from-xenograft) to separate the human cancer cell reads and mouse stroma cell reads.

We only use human specific reads or mouse specific reads for downstream analysis. One problem of **Xenome 1.0** is that the `+` symbols in **fastq** file are missing, you can use the command below to fix this problem and gzip the file simutaneously.

```
sed 's/^NS500589/@NS500589/' read.fastq | sed 's/^\s*$/+/' | gzip > read.fastq.gz
```

You can also use the custom python script `XenomeCorrect.py`.

```
XenomeCorrect.py read.fastq.fastq | gzip > read.fastq.gz
```

## Read Mapping

We used [STAR 2.4.1d](https://github.com/alexdobin/STAR) to map the RNA-seq reads. In order to improve mapping accuracy, the gene transfer format file was supplied at the genome index generation step with the command line option `--sjdbOverhang 79` (**ReadLength - 1**) together with the genome fasta file.

The bash code used for STAR reference index building and mapping can be found in **SupplementaryCode1.sh**.

## Reference Sample for CIBERSORT

We used [CIBERSORT](https://cibersort.stanford.edu/) to estimate the bone stroma percent in in bone-in-culture array (**BICA**) models and *in vivo* bone lesions (**IVBL**).

**CIBERSORT** only provide signature gene file for 22 immune cell types (LM22), and we need to generate a custom 13 signatures of bone marrow cell types using public available RNA-seq dataset [E-MTAB-2923](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2923/).

The directory **boneRefSamples** contains:

File Name  | Description
------------- | -------------
Ref folder | read quantification table for 13 samples
SupplementaryCode2.Rmd  | R code for generating signature gene file
BICA_DESeq.txt | signature gene file for 13 bone marrow cells
Samples\_Ortho\_BICA_Bone.txt | gene expression for RNA-seq samples (orthotopic, BICA, IVBL) in the project

You should be able to reproduce **Figure 1e** using the **BICA_DESeq.txt** and **Samples\_Ortho\_BICA_Bone.txt**.

## ssGSEA

**SupplementaryCode3.Rmd** is the code for generating ssGSEA input.

We run ssGSEA on [GenePattern](https://genepattern.broadinstitute.org/gp/pages/login.jsf). We use **ssGSEAProjection v7** module, not other later beta version.

We attached the outputs of ssGSEA in the folder **ssGSEAproj**.