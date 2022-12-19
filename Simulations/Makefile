# Adapted from: https://github.com/dobinlab/STARsoloManuscript.
# Generate simulated reads and TE counts from droplet-based (e.g. 10x Genomics)
# scRNA-seq data.
# Usage: `make all`

SHELL=/bin/bash
DATAURL := https://cg.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_fastqs.tar
RMSKURL :=  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz
GENOMEURL := http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
CBUMILEN := 26
CBLEN := 16
READLEN := 98
WLURL := https://github.com/10XGenomics/cellranger/raw/a83c753ce641db6409a59ad817328354fbe7187e/lib/python/cellranger/barcodes/737K-august-2016.txt
################################################################################

all: sims_counts.mtx sims_R1.fq sims_R2.fq

# prepare input for Seurat
sims_mex/: whitelist.txt genes.tsv sims_count.mtx
mkdir -p $@
cat whitelist.txt | gzip > $@/barcodes.tsv.gz
cat whitelist.txt | wc -l > barcodes_count.txt
cat genes.tsv | gzip > $@/features.tsv.gz
cat genes.tsv | wc -l > features_count.txt
awk '{s+=$$3} END {print s}' sims_count.mtx > umi_count.txt
echo '%%MatrixMarket matrix coordinate real general' | gzip > $@/matrix.mtx.gz
paste -d' ' barcodes_count.txt features_count.txt umi_count.txt | gzip >> $@/matrix.mtx.gz
cat sims_count.mtx | gzip >> $@/matrix.mtx.gz

# make count matrix
sims_counts.mtx: sims_R1.fq genes.tsv whitelist.txt
	awk -O -b -v cbLen=$(CBLEN) -v geneFile=genes.tsv \
		-v whiteList=whitelist.txt -f truthFromCBUMI.awk $< \
		| sort -k2,2n -k1,1n > $@

# link FASTQs to R1 R2 in samples/
sims_R2.fq: cDNA.fq
	ln -fs $< $@
sims_R1.fq: CBUMI.fq
	ln -fs $< $@

# generate FASTQs from bwa
cDNA.fq:
	$(MAKE) log.SIMS_FASTQ
CBUMI.fq:
	$(MAKE) log.SIMS_FASTQ
log.SIMS_FASTQ: genome_transcripts.1line.fa transcript_to_gene.txt genes.tsv whitelist.txt Mapped.txt
	set -v &&\
		awk -O -b -v errRate=0.005 -v readLen=$(READLEN) -v bcLen=$(CBUMILEN) \
			-v multiGene=yes -f CBUMIcDNAfqFromMapped_GenicGenomic.awk \
			$^ >& $@

whitelist.txt:
	wget -q -O $@ $(WLURL)

genes.tsv: transcript_to_gene.txt
	awk '{if (!($$2 in G)) print $$2; G[$$2]=0}' $< > $@

#transcript_to_gene.2col.txt: transcript_to_gene.txt
#	awk '{print $$1 "\t" $$2}' $< > $@

transcript_to_gene.txt: annotations.gtf
	awk '$$3=="exon" {gene=$$0; gsub(/.*gene_id "/,"",gene); \
		gsub(/".*/,"",gene); tr=$$0; gsub(/.*transcript_id "/,"",tr); \
		gsub(/".*/,"",tr); print tr "\t" gene }' $< > $@

# Extract mapped reads from SAM. Also unmapped, but they are not used.
Mapped.txt: aligned.bam
	samtools view -h $< |\
		awk -O -b -v bcLen=$(CBUMILEN) -v cbLen=$(CBLEN) \
		-f samToMappedAndUnmapped_mult.awk

# map sample reads to genome+transcriptome
aligned.bam: bwa_index/ cDNA_CBUMI.fq.gz
	bwa mem -a -M -t 20 $^ 2>$@.log | samtools view -b -o $@ -

# record CBUMI read from real data into readID of cDNA read for bwa
cDNA_CBUMI.fq.gz: samples/
	paste <(zcat $<fastqs/*_R1_*) <(zcat $<fastqs/*_R2_*) \
		| awk -O -b '{getline; cbumi=substr($$1,1,$(CBUMILEN));\
			cdna=$$2; getline;getline; ii++; name= "@" ii "_" cbumi;\
			print name "\n" cdna "\n+\n" $$2}' \
		| gzip > $@

# make bwa index
bwa_index/: genome_transcripts.fa
	mkdir $@ && bwa index -p $@ $< 2>$@.log

# concatenate genome and transcriptome
genome_transcripts.fa: genome.fa transcripts.fa
	cat $^ > $@

# concatenate genome and transcriptome on a single line
genome_transcripts.1line.fa: genome.fa transcripts.fa
	cat $^ \
		| awk -O -b '{if ($$1~/^>/) $$0="\n" substr($$1,2) " "; printf $$0}' > $@

# extract sequence from annotation
transcripts.fa: genome.fa genome.fa.fai annotations.gtf
	gffread -w $@ -g genome.fa annotations.gtf

# filter gtf by chromosome and replace "_" to "-" in gene or transcript names
annotations.gtf: rmsk.gtf genome.fa.fai
	zcat $< | awk -F"\t" -vOFS="\t" ' \
		FNR==NR { a[$$1];next } \
		$$1 in a { gsub("_","-",$$9); gsub("-id ","_id ",$$9); print }' genome.fa.fai - > $@

# download and index genome fasta
genome.fa.fai: genome.fa
	samtools faidx $<
genome.fa:
	curl -s $(GENOMEURL) | zcat > $@

# download repeatmasker annotation and convert to gtf
rmsk.gtf:
	wget -q $(RMSKURL) && zcat hg38.fa.out.gz | ./rmsk2bed -f gtf -t > $@

# download fastqs
samples/:
	mkdir $@ && cd $@ && wget -q $(DATAURL) && tar -xf *.tar
