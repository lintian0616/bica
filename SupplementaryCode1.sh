# Download Genome Reference
cd $GenomesRef
wget ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.83.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.83.chr.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-83/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-83/gtf/mus_musculus/Mus_musculus.GRCm38.83.chr.gtf.gz
gunzip Mus_musculus.GRCm38.83.chr.gtf.gz

# In silico Sorting Reads Sorting
cd $XenomeIdx
xenome-1.0.1-r/xenome index -M 8 -T 8 -P kmer25 --tmp-dir ./ -v -H Mus_musculus.GRCm38.dna.primary_assembly.fa -G Homo_sapiens.GRCh38.dna.primary_assembly.fa
xenome-1.0.1-r/xenome classify -P $XenomeIdx --pairs -i in_1.fastq -i in_2.fastq --graft-name human --host-name mouse

# Processing 2D/Orthotopic/BICA/Bone RNA-seq Samples
cd $RSEM_Idx
## For Batch 1: set sjdboverhang 79; for Batch 2: set sjdboverhang 75 (read length - 1)
rsem-prepare-reference -p 8 --gtf /$GenomesRef/Homo_sapiens.GRCh38.83.chr.gtf --star --star-sjdboverhang 79 $Genomesref/Homo_sapiens.GRCh38.primary.fa GRCh38
rsem-calculate-expression --star -p 8 --star-gzipped-read-file --output-genome-bam --calc-ci --estimate-rspd --quiet --phred33-quals --paired-end human_1.fastq.gz human_2.fastq.gz $RSEM_Idx/GRCh38 sample_human_RSEM

# Processing Bone Marrow Reference Samples
## Build Genome Index
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $GenomesDir/Mouse_STAR_69 --genomeFastaFiles $GenomesRef/Mus_musculus.GRCm38.dna.primary_assembly.fa --limitGenomeGenerateRAM 120000000000 --sjdbGTFfile $GenomeRef/Mus_musculus.GRCm38.83.chr.gtf --sjdbOverhang 69
## Align Reads to Using STAR
STAR --genomeDir $GenomesRef/Mouse_STAR_69 --readFilesIn sample_1.fastq.gz sample_2.fastq.gz --readFilesCommand zcat --outSAMstrandField intronMotif --runThreadN 8 --outFileNamePrefix sample_
## Sort the Aligned Reads
samtools view -bS -h sample_Aligned.out.sam | samtools sort -n -@ 8 -m 8G - sample_sorted
## Quantification Using HTseq
htseq-count -f bam -m union --strand=no --order=name sample_sorted.bam $GenomesRef/GRCm38.83.chr.gtf > sample_htseq.txt

