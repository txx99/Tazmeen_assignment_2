# Variables
SRA = "SRR1972739"
REF_ID = "AF086833.2"
RESULTS_FOLDER = "results"
RAW_DIR=f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR=f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR=f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR=f"{RESULTS_FOLDER}/annotated"
QC_DIR=f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR=f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR=f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR=f"{RESULTS_FOLDER}/snakemake"
BUCKET="sohail-binf5506"
S3_PREFIX="ebola"
 
rule all:
    input: 
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/reference.dict",
        f"{RAW_DIR}/reference.fasta.fai",
        f"{ALIGNED_DIR}/aligned.sam",
        f"{ALIGNED_DIR}/aligned.sorted.bam",
        f"{ALIGNED_DIR}/.validated_bam",
        f"{ALIGNED_DIR}/dedup.bam",
        f"{ALIGNED_DIR}/dup_metrics.txt",
        f"{ALIGNED_DIR}/dedup.bam.bai",
        f"{VARIANT_DIR}/raw_variants.vcf",
        f"{VARIANT_DIR}/filtered_variants.vcf",
        f"{SNPEFF_DATA_DIR}/genes.gbk",
        f"{SNPEFF_DIR}/snpEff.config",
        f"{SNPEFF_DATA_DIR}/snpEffectPredictor.bin"
        # f"{SNAKEMAKE_DIR}/.s3_upload_done"
 
rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        """
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """
 
rule download_reference:
    input:
        marker = rules.create_dirs.output.marker
    output:
        reference_fasta = f"{RAW_DIR}/reference.fasta"
    shell:
        """
        echo Downloading reference genome...
        efetch -db nucleotide -id {REF_ID} -format fasta > {RAW_DIR}/reference.fasta
        echo Downloaded reference genome!
        """

rule download_sra:
    input:
        marker = rules.create_dirs.output.marker
    output:
        sequence_sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        """
        echo Downloading sequencing data...
        prefetch {SRA} -O {RAW_DIR}
        echo Downloaded sequencing data!
        """
 
rule extract_sequence:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_sra = rules.download_sra.output.sequence_sra
    output:
        sequence_fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        echo Extracting sequencing data...
        fastq-dump -X 10000 {RAW_DIR}/{SRA}/{SRA}.sra -O {RAW_DIR}
        echo Extracted sequencing data!
        """

rule run_fastqc:
    input:
        marker = rules.create_dirs.output.marker,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq
    output:
        fastqc = f"{QC_DIR}"
    shell:
        """
        echo Running FastQC on raw reads...
        fastqc -o {QC_DIR} {RAW_DIR}/{SRA}.fastq 
        """
 
rule indexing:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta
    output:
        fai = f"{RAW_DIR}/reference.fasta.fai"
    shell:
        """
        echo Indexing reference genome with samtools...
        samtools faidx {RAW_DIR}/reference.fasta
        """


rule fasta_dict:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        fai = rules.indexing.output.fai
    output:
        fa_dict = f"{RAW_DIR}/reference.dict"
    shell:
        """
        echo Creating FASTA dictionary using GATK...
        gatk CreateSequenceDictionary -R {RAW_DIR}/reference.fasta -O {RAW_DIR}/reference.dict
        """

rule bwa_index:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta
    output:
        bwt = f"{RAW_DIR}/reference.fasta.bwt"  # Snakemake only needs one file to track
    shell:
        """
        echo Indexing reference genome with BWA...
        bwa index {input.reference_fasta}
        """


rule aligning:
    input: 
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq,
        bwt = rules.bwa_index.output.bwt,
        fai = rules.indexing.output.fai
    output:
        aligned_sam = f"{ALIGNED_DIR}/aligned.sam"
    shell:
        """
        echo Aligning reads with read groups...
        bwa mem -R '@RG\\tID:1\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:sample1' {input.reference_fasta} {input.sequence_fastq} > {output.aligned_sam}
        echo Aligned reads!
        """
    
rule bam_conversion:
    input:
        marker = rules.create_dirs.output.marker,
        aligned_sam = rules.aligning.output.aligned_sam
    output:
        bamed = f"{ALIGNED_DIR}/aligned.sorted.bam"
    shell:
        """
        echo Converting SAM to sorted BAM...
        samtools view -b {ALIGNED_DIR}/aligned.sam | samtools sort -o {ALIGNED_DIR}/aligned.sorted.bam
        """

rule bam_validaxn:
    input:
        marker = rules.create_dirs.output.marker,
        bamed = rules.bam_conversion.output.bamed
    output:
        valid_bam = f"{ALIGNED_DIR}/.validated_bam"
    shell:
        """
        echo Validating BAM file...
        gatk ValidateSamFile -I {ALIGNED_DIR}/aligned.sorted.bam -MODE SUMMARY > {ALIGNED_DIR}/.validated_bam
        """

rule mark_dups:
    input:
        marker = rules.create_dirs.output.marker,
        bamed = rules.bam_conversion.output.bamed
    output:
        dedup_bam = f"{ALIGNED_DIR}/dedup.bam",
        dup_mets = f"{ALIGNED_DIR}/dup_metrics.txt"
    shell:
        """
        echo Marking duplicates...
        gatk MarkDuplicates -I {ALIGNED_DIR}/aligned.sorted.bam -O {ALIGNED_DIR}/dedup.bam -M {ALIGNED_DIR}/dup_metrics.txt
        """

rule dedup_indexing:
    input:
        marker = rules.create_dirs.output.marker,
        dedup_bam = rules.mark_dups.output.dedup_bam
    output:
        dedup_bam_index = f"{ALIGNED_DIR}/dedup.bam.bai"
    shell:
        """
        echo Indexing deduplicated BAM file...
        samtools index {ALIGNED_DIR}/dedup.bam
        """

rule variant_cal:
    input: 
        marker = rules.create_dirs.output.marker,
        dedup_bam = rules.mark_dups.output.dedup_bam,
        reference_fasta = rules.download_reference.output.reference_fasta,
        dedup_bam_index = rules.dedup_indexing.output.dedup_bam_index
    output:
        raw_vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    shell:
        """
        echo Calling variants...
        gatk HaplotypeCaller -R {RAW_DIR}/reference.fasta -I {ALIGNED_DIR}/dedup.bam -O {VARIANT_DIR}/raw_variants.vcf
        echo Called variants!
        """

rule filter_variants:
    input: 
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        raw_variants = rules.variant_cal.output.raw_vcf
    output:
        filtered_variants = f"{VARIANT_DIR}/filtered_variants.vcf"
    shell:
        """
        echo Filtering variants...
        gatk VariantFiltration -R {RAW_DIR}/reference.fasta -V {VARIANT_DIR}/raw_variants.vcf -O {VARIANT_DIR}/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
        """

rule snpEff_ref_download:
    input:
        marker = rules.create_dirs.output.marker
    output: 
        genbank_ref = f"{SNPEFF_DATA_DIR}/genes.gbk"
    shell:
        """
        echo Downloading reference GenBank file for snpEff...
        efetch -db nucleotide -id {REF_ID} -format genbank > {SNPEFF_DATA_DIR}/genes.gbk
        echo Downloaded GenBank file for snpEff!
        """

rule snpEff_config:
    input:
        marker = rules.create_dirs.output.marker,
        reference_fasta = rules.download_reference.output.reference_fasta,
        genbank_ref = rules.snpEff_ref_download.output.genbank_ref
    output:
        snpEff_config_file = f"{SNPEFF_DIR}/snpEff.config"
    shell:
        """
        echo Creating custom snpEff configuration file...
        cat <<EOF > {SNPEFF_DIR}/snpEff.config
        # Custom snpEff config for reference_db
        reference_db.genome : reference_db
        reference_db.fa : $(readlink -f {RAW_DIR}/reference.fasta)
        reference_db.genbank : $(readlink -f {SNPEFF_DATA_DIR}/genes.gbk)
        EOF
        """

rule build_snpEff_db: # builds a folder --> selecting one file as output
    input:
        marker = rules.create_dirs.output.marker,
        snpEff_config_file = rules.snpEff_config.output.snpEff_config_file 
    output:
        snpEff_db = f"{SNPEFF_DATA_DIR}/snpEffectPredictor.bin"
    shell:
        """
        echo Building snpEff database...
        snpEff build -c {SNPEFF_DIR}/snpEff.config -genbank -v -noCheckProtein reference_db
        echo Built snpEff database!
        """



# # rule upload_s3:
# #     input:
# #         reference_fasta = rules.download_reference.output.reference_fasta,
# #         sequence_sra = rules.download_sra.output.sequence_sra,
# #         sequence_fastq = rules.extract_sequence.output.sequence_fastq
# #     output:
# #         marker = f"{SNAKEMAKE_DIR}/.s3_upload_done"
# #     run:
# #         import os
# #         import boto3
# #         s3 = boto3.client("s3")
 
# #         for root, dirs, files in os.walk(RESULTS_FOLDER):
# #             for file in files:
# #                 local_file = os.path.join(root, file)
# #                 relative_path = os.path.relpath(local_file, RESULTS_FOLDER)
# #                 s3_key = os.path.join(S3_PREFIX, relative_path).replace("\\", "/")
 
# #                 print(f"Uploading {local_file} to s3://{BUCKET}/{s3_key}")
# #                 s3.upload_file(local_file, BUCKET, s3_key)
 
# #         with open(output.marker, "w") as f:
# #             f.write("Upload Complete!")