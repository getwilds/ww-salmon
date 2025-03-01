version 1.0

workflow SalmonRnaSeq {
    input {
        # Essential inputs for paired-end data only
        File transcriptome_fasta  # Reference transcriptome in FASTA format
        Array[File] fastq_r1_files  # FASTQ files for read 1
        Array[File] fastq_r2_files  # FASTQ files for read 2
        
        # Docker image
        String salmon_docker = "combinelab/salmon:latest"
    }
    
    # Build Salmon index
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
                fastq_r2 = fastq_r2_files[i],
                docker_image = salmon_docker
        }
    }
    
    # Merge Salmon quant results
    call MergeSalmonResults {
        input:
            salmon_quant_dirs = SalmonQuant.salmon_quant_dir,
            sample_names = SalmonQuant.output_sample_name,
            docker_image = salmon_docker
    }
    
    output {
        File salmon_index_tar = BuildSalmonIndex.salmon_index
        Array[File] salmon_quant_dirs = SalmonQuant.salmon_quant_dir
        File merged_tpm_matrix = MergeSalmonResults.tpm_matrix
        File merged_counts_matrix = MergeSalmonResults.counts_matrix
        File sample_list = MergeSalmonResults.sample_list
    }
    
    parameter_meta {
        transcriptome_fasta: "Reference transcriptome in FASTA format"
        fastq_r1_files: "FASTQ files for read 1"
        fastq_r2_files: "FASTQ files for read 2"
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
        File fastq_r2
        String docker_image
        Int memory_gb = 16
        Int cpu = 8
        Int disk_size_gb = 100
    }
    
    command <<<
        set -e
        
        # Extract the Salmon index
        mkdir -p salmon_index
        tar -xzf ~{salmon_index_dir} -C ./
        
        # Paired-end quantification with best practice parameters
        salmon quant \
            -i salmon_index \
            --libType A \
            -1 ~{fastq_r1} \
            -2 ~{fastq_r2} \
            -o ~{sample_name}_quant \
            -p ~{cpu} \
            --validateMappings \
            --gcBias \
            --seqBias
        
        # Tar up the output directory
        tar -czf ~{sample_name}_quant.tar.gz ~{sample_name}_quant
    >>>
    
    output {
        File salmon_quant_dir = "~{sample_name}_quant.tar.gz"
        String output_sample_name = "~{sample_name}"
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
            base_name=$(basename $quant_dir .tar.gz)
            mkdir -p "quant_dirs/$base_name"
            tar -xzf $quant_dir -C quant_dirs/
        done
        
        # Create a list of sample names for reference
        echo "~{sep="\n" sample_names}" > sample_names.txt
        
        # Create a list of quant directories for quantmerge
        quant_dirs_list=()
        for sample in ~{sep=" " sample_names}; do
            quant_dirs_list+=("quant_dirs/${sample}_quant")
        done
        
        # Merge transcript-level TPM values
        salmon quantmerge \
            --quants "${quant_dirs_list[@]}" \
            --column TPM \
            -o tpm_matrix.tsv
        
        # Merge transcript-level estimated count values
        salmon quantmerge \
            --quants "${quant_dirs_list[@]}" \
            --column NumReads \
            -o counts_matrix.tsv
    >>>
    
    output {
        File tpm_matrix = "tpm_matrix.tsv"
        File counts_matrix = "counts_matrix.tsv"
        File sample_list = "sample_names.txt"
    }
    
    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: 1
    }
}
