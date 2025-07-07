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
        f"{SNAKEMAKE_DIR}/.s3_upload_done"
 
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
    output:
        index = {RAW_DIR}/reference.fasta
    shell:
        """
        echo Indexing reference genome with samtools...
        samtools faidx {RAW_DIR}/reference.fasta
        """







rule upload_s3:
    input:
        reference_fasta = rules.download_reference.output.reference_fasta,
        sequence_sra = rules.download_sra.output.sequence_sra,
        sequence_fastq = rules.extract_sequence.output.sequence_fastq
    output:
        marker = f"{SNAKEMAKE_DIR}/.s3_upload_done"
    run:
        import os
        import boto3
        s3 = boto3.client("s3")
 
        for root, dirs, files in os.walk(RESULTS_FOLDER):
            for file in files:
                local_file = os.path.join(root, file)
                relative_path = os.path.relpath(local_file, RESULTS_FOLDER)
                s3_key = os.path.join(S3_PREFIX, relative_path).replace("\\", "/")
 
                print(f"Uploading {local_file} to s3://{BUCKET}/{s3_key}")
                s3.upload_file(local_file, BUCKET, s3_key)
 
        with open(output.marker, "w") as f:
            f.write("Upload Complete!")