#!/usr/bin/env nextflow

params.star            = params.genomes[ params.genome ]?.star
params.repeatmasker_te = '../annotations.bed'
params.whitelist   = '../whitelist.txt'

ch_rmsk      = file(params.repeatmasker_te, checkIfExists: true)
ch_whitelist = file(params.whitelist, checkIfExists: true)

ch_reads = Channel
  .fromFilePairs('../sims_R{1,2}.fq', flat: true)
  .map { ['sim_pbmc8k', file(it[1]), file(it[2])] }


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

process STAR {
  cpus 6
  memory 50.GB
  publishDir "./STAR", mode: 'copy'

  input:
  tuple val(sample), path(read1), path(read2) from ch_reads
  path index from ch_star_index
  path gtf from ch_gtf
  path whitelist from ch_whitelist

  output:
  tuple val(sample), path('*.bam'), path('*Solo.out/Gene/filtered/barcodes.tsv') into ch_alignments, ch_alignments_scte
  tuple val(sample), path('*.bam') into ch_alignments_samtools
  path '*.{bam,out,tab}'

  script:
  def r1 = read1 instanceof List ? read1.join(',') : read1
  def r2 = read2 instanceof List ? read2.join(',') : read2
  def zcat = r1
  """
  STAR \\
    --runThreadN $task.cpus \\
    --genomeDir $index \\
    --outFileNamePrefix ${sample}. \\
    --sjdbGTFfile $gtf \\
    --readFilesIn $r2 $r1 \\
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

process SAMTOOLS_INDEX {
  cpus 1
  memory 1.GB
  publishDir "./STAR", mode: 'copy'

  input:
  tuple val(sample), path(bam) from ch_alignments_samtools

  output:
  path '*.bai' into ch_samtools_index

  script:
  """
  samtools index $bam
  """
}

ch_alignments_irescue = ch_alignments.join(ch_samtools_index)

process IRESCUE {
  cpus 8
  memory 16.GB
  publishDir "./IRESCUE", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(whitelist), path(bai) from ch_alignments_irescue

  output:
  path '*'

  script:
  """
  irescue -b $bam -g hg38 -w $whitelist -p $task.cpus \
    --CBtag CB --UMItag UR --keeptmp
  """
}

process SCTE {
  cpus 6
  memory 160.GB
  publishDir "./scTE", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(whitelist) from ch_alignments_scte
  path rmsk from ch_rmsk
  path gtf from ch_gtf

  output:
  path '*'

  script:
  """
  samtools view -D CB:$whitelist -o ${sample}_filtered.bam -@ 4 $bam
  scTE_build -te $rmsk -gene $gtf -o scte_index
  scTE -i ${sample}_filtered.bam -x scte_index -o $sample \\
    -CB CB -UMI UR --hdf5 True -p $task.cpus
  """
}
