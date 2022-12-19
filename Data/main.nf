#!/usr/bin/env nextflow

// Download GTF and convert chromosome names to UCSC format
process GET_GTF {
  cpus 1
  memory 1.GB

  output:
  path 'annotation.ucsc_names.gtf' into ch_gtf

  script:
  """
  curl --output annotation.gtf.gz --silent $params.gtf_url
  curl --silent $params.gencode_assembly_report \\
    | tr '\\015' \\
    | grep -v '^#' > report.txt
  zcat annotation.gtf.gz \\
    | gawk -F"\\t" '\\
      FNR==NR { chr[\$5]=\$10; next } \\
      { if(\$1 in chr) \$1=chr[\$1] } \\
      { OFS="\\t"; print }' > annotation.ucsc_names.gtf
  """
}

// Download genome fasta file
process GET_GENOME {
  cpus 1
  memory 1.GB

  output:
  path 'hg38.fa' into ch_genome

  script:
  """
  curl --output hg38.fa.gz --silent $params.genome_url
  gzip -d hg38.fa.gz
  """
}

// Build STAR index
process STAR_INDEX {
  cpus 4
  memory 36.GB

  input:
  path genome from ch_genome
  path gtf from ch_gtf

  output:
  path 'STARINDEX' into ch_star_index

  script:
  """
  STAR \\
    --runMode genomeGenerate \\
    --genomeDir STARINDEX \\
    --genomeFastaFiles $genome \\
    --sjdbGTFfile $gtf
  """
}

// Get URLs of fastq files from ArrayExpress samlple sheet
process GET_READS_URLS {
  cpus 1
  memory 100.MB

  output:
  path 'urls.csv' into ch_reads_urls

  script:
  """
  curl --output sdrf.txt --silent $params.sdrf_url
  gawk -F"\\t" 'NR>1 && \$11~/KUL/ {
    gsub(/.+\\(|\\)/,"",\$11)
    if (\$NF~/core/) {
      type="T"
    } else if (\$NF~/border/) {
      type="B"
    } else {
      type="N"
    }
    name=\$11"-"type
    OFS=","; print name,\$61,\$63
  }' sdrf.txt > urls.csv
  """
}

// Download fastqs
process GET_READS {
  cpus 1
  memory 100.MB

  input:
  val(sample), val(urlR1), val(urlR2) from ch_reads_urls.splitCsv()

  output:
  val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz") into ch_reads

  script:
  """
  curl --output ${sample}_R1.fastq.gz --silent $urlR1
  curl --output ${sample}_R2.fastq.gz --silent $urlR2
  """
}

// Download 10x Genomics whitelist
process GET_WHITELIST {
  cpus 1
  memory 100.MB

  output:
  path 'whitelist.txt' into ch_whitelist

  script:
  """
  curl --output whitelist.txt --silent $params.whitelist_url
  """
}

// Run STARSolo
process STAR {
  cpus 6
  memory 36.GB
  publishDir "${params.outdir}/${sample}/STAR", mode: 'copy'

  input:
  tuple val(sample), path(read1), path(read2) from ch_reads.groupTuple()
  path index from ch_star_index
  path whitelist from ch_whitelist

  output:
  tuple val(sample), path('*.bam'), path('*Solo.out/Gene/filtered/barcodes.tsv') into ch_alignments
  path '*.{bam,out,tab}'

  script:
  def r1 = read1 instanceof List ? read1.join(',') : read1
  def r2 = read2 instanceof List ? read2.join(',') : read2
  """
  STAR \\
    --runThreadN $task.cpus \\
    --genomeDir $index \\
    --outFileNamePrefix ${sample}. \\
    --readFilesIn $r2 $r1 \\
    --readFilesCommand zcat \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMattributes NH HI AS nM NM MD jM jI XS MC ch cN CR CY UR UY GX GN CB UB sM sS sQ \\
    --outFilterMultimapNmax 100 \\
    --winAnchorMultimapNmax 100 \\
    --twopassMode Basic \\
    --soloType CB_UMI_Simple \\
    --soloCellFilter EmptyDrops_CR 5000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
    --soloCBwhitelist $whitelist \\
    --limitOutSJcollapsed 2000000
  """
}

// Run IRescue
process IRESCUE {
  cpus 8
  memory 16.GB
  publishDir "${params.outdir}/${sample}", mode: 'copy', saveAs: { dir -> 'IRESCUE' }

  input:
  tuple val(sample), path(bam), path(whitelist) from ch_alignments
  path rmsk from ch_rmsk

  output:
  path 'IRescue_out'

  script:
  """
  samtools index $bam
  irescue -b $bam -g hg38 -w $whitelist -p $task.cpus --CBtag CB --UMItag UR
  """
}

