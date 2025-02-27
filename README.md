# ww-salmon RNA-seq Workflow
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)

## Overview

This workflow performs RNA-seq quantification using [Salmon](https://combine-lab.github.io/salmon/), a fast and accurate tool for transcript expression estimation. The workflow is designed to be simple to use while implementing best practices for RNA-seq analysis.

## What this Workflow Does

1. **Builds a Salmon Index** from your reference transcriptome
2. **Quantifies Transcripts** for each of your RNA-seq samples
3. **Generates Expression Matrices** combining results from all samples

## Requirements

- [Cromwell](https://github.com/broadinstitute/cromwell) or another WDL-compatible workflow engine
- Docker (the workflow uses the `combinelab/salmon` container)
- Input files:
  - Reference transcriptome (FASTA format)
  - RNA-seq reads (FASTQ format, can be gzipped)

## Quick Start

1. **Download the WDL file from this repository**:

2. **Create an inputs JSON file** (e.g., `inputs.json`):
   ```json
   {
     "SalmonRnaSeq.transcriptome_fasta": "path/to/transcriptome.fa",
     "SalmonRnaSeq.fastq_r1_files": [
       "path/to/sample1_R1.fastq.gz",
       "path/to/sample2_R1.fastq.gz"
     ],
     "SalmonRnaSeq.fastq_r2_files": [
       "path/to/sample1_R2.fastq.gz",
       "path/to/sample2_R2.fastq.gz"
     ]
   }
   ```

3. **Run the workflow** with Cromwell:
   ```
   java -jar cromwell.jar run salmon_rnaseq.wdl -i inputs.json
   ```

## Input Parameters

| Parameter | Description | Required? |
|-----------|-------------|-----------|
| `transcriptome_fasta` | Reference transcriptome in FASTA format | Yes |
| `fastq_r1_files` | Array of FASTQ files for read 1 (or single-end reads) | Yes |
| `fastq_r2_files` | Array of FASTQ files for read 2 (for paired-end data) | No |
| `salmon_docker` | Docker image for Salmon (default: "combinelab/salmon:latest") | No |

## Outputs

| Output | Description |
|--------|-------------|
| `salmon_index_tar` | Compressed Salmon index (can be reused for future analyses) |
| `salmon_quant_dirs` | Compressed quantification results for each sample |
| `merged_tpm_matrix` | Combined TPM values matrix for all samples |
| `merged_counts_matrix` | Combined read counts matrix for all samples |

## Examples

### For Paired-End Data

```json
{
  "SalmonRnaSeq.transcriptome_fasta": "references/gencode.v38.transcripts.fa",
  "SalmonRnaSeq.fastq_r1_files": [
    "samples/sample1_R1.fastq.gz",
    "samples/sample2_R1.fastq.gz"
  ],
  "SalmonRnaSeq.fastq_r2_files": [
    "samples/sample1_R2.fastq.gz", 
    "samples/sample2_R2.fastq.gz"
  ]
}
```

### For Single-End Data

```json
{
  "SalmonRnaSeq.transcriptome_fasta": "references/gencode.v38.transcripts.fa",
  "SalmonRnaSeq.fastq_r1_files": [
    "samples/sample1.fastq.gz",
    "samples/sample2.fastq.gz"
  ]
}
```

## Common Questions

### How do I get a transcriptome file?

You can download reference transcriptomes from:
- [GENCODE](https://www.gencodegenes.org/human/) (human/mouse)
- [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) (many species)
- [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables)

For human, a common choice is the GENCODE reference:
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
gunzip gencode.v38.transcripts.fa.gz
```

### How do I extract the results?

The workflow provides compressed output directories for each sample. To extract a specific sample's results:

```bash
tar -xzf sample1_quant.tar.gz
```

This will create a directory `sample1_quant` containing Salmon's output files, including:
- `quant.sf`: The main quantification results file
- `logs/`: Directory containing Salmon log files
- `lib_format_counts.json`: Information about the library type

### What are TPM and counts?

- **TPM (Transcripts Per Million)**: Normalized expression values suitable for comparing expression levels between samples
- **counts**: Estimated number of fragments/reads from each transcript, suitable for differential expression analysis

## Under the Hood

This workflow:
1. Creates a Salmon index from your transcriptome
2. Processes each sample with optimal settings:
   - Automatic library type detection
   - GC bias correction
   - Sequence-specific bias correction
   - Mapping validation
3. Combines results into unified matrices with properly labeled sample names

## Troubleshooting

**Error**: "Docker image not found"
- Solution: Ensure Docker is installed and running

**Error**: "File not found"
- Solution: Check the paths in your inputs.json file

**Error**: "Memory allocation failed"
- Solution: Adjust the `memory_gb` parameters in the WDL file

## Advanced Customization

If you need to modify the workflow for advanced settings:

1. Edit the runtime parameters at the task level:
   ```wdl
   runtime {
       docker: docker_image
       memory: "~{memory_gb} GB"
       cpu: cpu
       disks: "local-disk ~{disk_size_gb} SSD"
       preemptible: 1
   }
   ```

2. Add additional Salmon parameters in the command sections if needed

## Need Additional Help?

- [Salmon Documentation](https://salmon.readthedocs.io/en/latest/)
- [WDL Documentation](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)
- [Cromwell Documentation](https://cromwell.readthedocs.io/en/stable/)
