version 1.0

workflow SalmonRnaSeq {
    input {
        # Essential inputs only
        File transcriptome_fasta  # Reference transcriptome in FASTA format
        Array[File] fastq_r1_files  # FASTQ files for read 1 (or single-end reads)
        Array[File]? fastq_r2_files  # Optional FASTQ files for read 2 (paired-end)
        
        # Docker image
        String salmon_docker = "combinelab/salmon:latest"
    }
    
    # Check if paired-end or single-end data
    Boolean is_paired_end = defined(fastq_r2_files) && length(select_first([fastq_r2_files, []])) > 0
    
    # Build Salmon index (with simplified parameters)
    call BuildSalmonIndex {
        input:
            transcriptome_fasta = transcriptome_fasta,
            docker_image = salmon_docker
    }
    
    # Quantify each sample with Salmon
    scatter (i in range(length(fastq_r1_files))) {
        call SalmonQuant {
            input:
                salmon_index_dir = BuildSalmonIndex.salmon_index,
                sample_name = basename(fastq_r1_files[i], ".fastq.gz"),
                fastq_r1 = fastq_r1_files[i],
                fastq_r2 = if is_paired_end then select_first([fastq_r2_files, []])[i] else null,
                docker_image = salmon_docker
        }
    }
    
    # Merge Salmon quant results
    call MergeSalmonResults {
        input:
            salmon_quant_dirs = SalmonQuant.salmon_quant_dir,
            sample_names = SalmonQuant.sample_name,
            docker_image = salmon_docker
    }
    
    output {
        # Simplified outputs
        File salmon_index_tar = BuildSalmonIndex.salmon_index
        Array[File] salmon_quant_dirs = SalmonQuant.salmon_quant_dir
        File merged_tpm_matrix = MergeSalmonResults.tpm_matrix
        File merged_counts_matrix = MergeSalmonResults.counts_matrix
    }
    
    parameter_meta {
        transcriptome_fasta: "Reference transcriptome in FASTA format"
        fastq_r1_files: "FASTQ files for read 1 (or single-end reads)"
        fastq_r2_files: "Optional FASTQ files for read 2 (paired-end)"
        salmon_docker: "Docker image for Salmon"
    }
}

task BuildSalmonIndex {
    input {
        File transcriptome_fasta
        String docker_image
        Int memory_gb = 16
        Int cpu = 8
        Int disk_size_gb = 100
    }
    
    command <<<
        set -e
        
        # Build Salmon index with default parameters and decoys
        salmon index \
            -t ~{transcriptome_fasta} \
            --tmpdir=./tmp \
            -i salmon_index \
            -p ~{cpu} \
            --gencode
        
        # Create tar archive of the index
        tar -czf salmon_index.tar.gz salmon_index
    >>>
    
    output {
        File salmon_index = "salmon_index.tar.gz"
    }
    
    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: 1
    }
}

task SalmonQuant {
    input {
        File salmon_index_dir
        String sample_name
        File fastq_r1
        File? fastq_r2
        String docker_image
        Int memory_gb = 16
        Int cpu = 8
        Int disk_size_gb = 100
    }
    
    Boolean is_paired_end = defined(fastq_r2)
    
    command <<<
        set -e
        
        # Extract the Salmon index
        mkdir -p salmon_index
        tar -xzf ~{salmon_index_dir} -C ./
        
        # Paired-end or single-end quantification with common best practice parameters
        if [ "~{is_paired_end}" == "true" ]; then
            # Paired-end quantification
            salmon quant \
                -i salmon_index \
                -1 ~{fastq_r1} \
                -2 ~{fastq_r2} \
                --libType A \
                -o ~{sample_name}_quant \
                -p ~{cpu} \
                --validateMappings \
                --gcBias \
                --seqBias
        else
            # Single-end quantification
            salmon quant \
                -i salmon_index \
                -r ~{fastq_r1} \
                --libType A \
                -o ~{sample_name}_quant \
                -p ~{cpu} \
                --validateMappings \
                --gcBias \
                --seqBias
        fi
        
        # Tar up the output directory
        tar -czf ~{sample_name}_quant.tar.gz ~{sample_name}_quant
    >>>
    
    output {
        File salmon_quant_dir = "~{sample_name}_quant.tar.gz"
        String sample_name = "~{sample_name}"
    }
    
    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: 1
    }
}

task MergeSalmonResults {
    input {
        Array[File] salmon_quant_dirs
        Array[String] sample_names
        String docker_image
        Int memory_gb = 8
        Int cpu = 4
        Int disk_size_gb = 50
    }
    
    command <<<
        set -e
        
        # Extract all quantification results
        mkdir -p quant_dirs
        for quant_dir in ~{sep=" " salmon_quant_dirs}; do
            tar -xzf $quant_dir -C quant_dirs/
        done
        
        # Merge transcript-level TPM values
        salmon quantmerge \
            --quants quant_dirs/*/quant.sf \
            --column TPM \
            --genes \
            -o tpm_matrix.tsv
        
        # Merge transcript-level estimated count values
        salmon quantmerge \
            --quants quant_dirs/*/quant.sf \
            --column NumReads \
            --genes \
            -o counts_matrix.tsv
        
        # Create a list of sample names for column headers
        echo "~{sep="\n" sample_names}" > sample_names.txt
        
        # Rename columns with sample names for better readability
        python3 - <<EOF
import pandas as pd
import os

# Read the current matrices
tpm = pd.read_csv("tpm_matrix.tsv", sep="\t")
counts = pd.read_csv("counts_matrix.tsv", sep="\t")

# Read sample names
sample_names = [line.strip() for line in open("sample_names.txt").readlines()]

# Get the current column names (excluding 'Name')
tpm_cols = list(tpm.columns)[1:]
counts_cols = list(counts.columns)[1:]

# Ensure we have the right number of samples
if len(tpm_cols) == len(sample_names):
    # Create a mapping of old column names to sample names
    tpm_colmap = {old: new for old, new in zip(tpm_cols, sample_names)}
    counts_colmap = {old: new for old, new in zip(counts_cols, sample_names)}
    
    # Rename the columns
    tpm = tpm.rename(columns=tpm_colmap)
    counts = counts.rename(columns=counts_colmap)
    
    # Save the updated matrices
    tpm.to_csv("tpm_matrix_renamed.tsv", sep="\t", index=False)
    counts.to_csv("counts_matrix_renamed.tsv", sep="\t", index=False)
    
    # Move renamed files to final output
    os.rename("tpm_matrix_renamed.tsv", "tpm_matrix.tsv")
    os.rename("counts_matrix_renamed.tsv", "counts_matrix.tsv")
else:
    print(f"Warning: Number of samples ({len(sample_names)}) doesn't match number of columns ({len(tpm_cols)})")
EOF
    >>>
    
    output {
        File tpm_matrix = "tpm_matrix.tsv"
        File counts_matrix = "counts_matrix.tsv"
    }
    
    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: 1
    }
}
