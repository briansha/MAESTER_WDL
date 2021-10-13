version 1.0

## Version 10-13-2021
##
## This workflow runs MAESTER.
## MAESTER Documentation: https://github.com/vangalenlab/MAESTER-2021
##
## Note: Maegatk task requires bcall or support parameter to be set to "true".
##
## Cromwell version support - Successfully tested on v68
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow maester {
    input {
        String docker = "briansha/maester:4.1.0"                  # R 4.1.0, Ubuntu 18.04, and some R packages. - Used for all tasks using R scripts.
        String docker_homer = "briansha/maester_homer:4.11"       # v4.11 of HOMER
        String docker_star = "briansha/star:2.7.9"                # v2.7.9 of STAR
        String docker_samtools = "briansha/maester_samtools:1.13" # v1.13 of samtools.
        String docker_maegatk = "briansha/maester_maegatk:v01"    # v01 of maegatk.
    }

    call FilterBarcodes {
      input:
        docker = docker
    }

    call TrimWithHomer {
      input:
        docker = docker_homer,
        fastq_read2 = FilterBarcodes.output_fastq
    }

    call Star {
      input:
        docker = docker_star,
        fastq_read2 = TrimWithHomer.output_fastq
    }

    call TagCbUmi {
      input:
        docker = docker_samtools,
        input_sam = Star.output_sam
    }

    call SubsetForChrM {
      input:
        docker = docker_samtools,
        input_bam_10x = TagCbUmi.output_bam_10x
    }

    call MergeBamFiles {
      input:
        docker = docker_samtools,
        input_maester_bam = SubsetForChrM.output_bam,
        input_maester_bai = SubsetForChrM.output_bai
    }

    call Maegatk {
      input:
        docker = docker_maegatk,
        input_bam = MergeBamFiles.output_merged_bam,
        input_bai = MergeBamFiles.output_merged_bai
    }

    call MtCoverage {
      input:
        docker = docker,
        input_folder = Maegatk.output_zipped_file
    }

    output {
        File maegatk_output = Maegatk.output_zipped_file
        File mt_coverage_plots = MtCoverage.output_plot
    }

    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "This workflow runs MAESTER: https://github.com/vangalenlab/MAESTER-2021"
    }
}

# Filter for cell barcodes (CBs) and generate fastq files with CB
# and unique molecular identifiers (UMIs) from the Read 1 fastq in the read ID of the Read 2 fastq
task FilterBarcodes {
    input {
        File r_script        # Assemble_fastq.R - from the MAESTER github.

        # Maester 1.1 parameters
        File folder          # One or multiple directories containing fastq files, not searched recursively
        String sampleName    # Sample name that will be used for output files
        File cellBarcodes    # Allowlist of cell barcodes to filter by. Could be cells from the CellRanger filtered_feature_bc_matrix, or that passed scRNA-seq QC, or all whitelisted 10x cell barcodes.
        Int CBlength         # Length of the cell barcode (16 for 10x 3' v3)
        Int UMIlength        # Length of the UMI (12 for 10x 3' v3)

        # Runtime
        String docker
        Int? disk_size_override
        Int cpu = 1
        Float memory = 32.0
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float folder_size = size(folder, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * folder_size)])
    String folderName = basename(folder, ".zip")

    # Useful commands
    # ls
    # du -d 1 -h
    command <<<
        set -euo pipefail
        unzip ~{folder}

        R < ~{r_script} --no-save --args ~{folderName} ~{sampleName} ~{cellBarcodes} ~{CBlength} ~{UMIlength}
    >>>

    output {
        File output_fastq = "${sampleName}.fastq.gz"
        File output_stats = "${sampleName}.stats.txt"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Trim 24 bp from the start of Read 2
# (uninformative primer binding site)
task TrimWithHomer {
    input {
        File fastq_read2

        # Runtime
        String docker
        Int? disk_size_override
        Int cpu = 1
        Float memory = 3.5
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float fastq_read2_size = size(fastq_read2, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 4.0 * fastq_read2_size)])
    String fastq_read2_filename = basename(fastq_read2)

    command <<<
        set -euo pipefail
        mv ~{fastq_read2} .

        homerTools trim -5 24 ~{fastq_read2_filename}
    >>>

    output {
        File output_fastq = "${fastq_read2_filename}"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Align to hg38 reference genome using STAR.
task Star {
    input {
        File fastq_read2                       # Need to unzip this file if it is a .gz file for STAR to work.
        File genomeDir                         # Zip file for hg38 genome indexes from STAR
        String output_name = "Aligned.out.sam"

        # Runtime
        String docker
        Int? disk_size_override
        Int? cpu_override
        Float memory = 32.0
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float genomeDir_size = size(genomeDir, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 7.0 * genomeDir_size)])
    Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * floor(memory / 8) else 1])
    String fastq_read2_filename = basename(fastq_read2)
    String fastq_read2_unzipped_filename = basename(fastq_read2, ".gz")

    command <<<
        set -euo pipefail
        unzip ~{genomeDir}
        mv ~{fastq_read2} .
        gunzip ~{fastq_read2_filename}

        STAR \
        --runThreadN ~{cpu} \
        --genomeDir star \
        --readFilesIn ~{fastq_read2_unzipped_filename}
    >>>

    output {
        File output_sam = "${output_name}"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Add CB and UMI as sam tags.
task TagCbUmi {
    input {
        File input_sam
        File bash_script  # Tag_CB_UMI.sh - from the MAESTER github.

        # Runtime
        String docker
        Int? disk_size_override
        Int cpu = 1
        Float memory = 3.5
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float input_sam_size = size(input_sam, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 3.0 * input_sam_size)])
    String input_bam = basename(input_sam, ".sam") + ".bam"
    String output_bam_name = basename(input_sam, ".sam") + ".10x.bam"
    String bash_script_current_dir = basename(bash_script)

    #Convert sam to bam - then execute Tag_CB_UMI.sh
    command <<<
        set -euo pipefail

        samtools view -S -b ~{input_sam} > ~{input_bam}
        mv ~{bash_script} .
        chmod 700 ~{bash_script_current_dir}
        ./~{bash_script_current_dir} ~{input_bam}
    >>>

    output {
        File output_bam_10x = "${output_bam_name}"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Subset for chrM.
task SubsetForChrM {
    input {
        File input_bam_10x

        # Runtime
        String docker
        Int? disk_size_override
        Int cpu = 1
        Float memory = 3.5
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float input_bam_10x_size = size(input_bam_10x, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 3.0 * input_bam_10x_size)])
    String input_bam_10x_sorted = basename(input_bam_10x, ".bam") + ".sorted.bam"

    command <<<
        set -euo pipefail

        samtools sort ~{input_bam_10x} -o ~{input_bam_10x_sorted}
        samtools index ~{input_bam_10x_sorted} ~{input_bam_10x_sorted}.bai
        samtools view ~{input_bam_10x_sorted} chrM -b > ~{input_bam_10x_sorted}.subset.bam
        samtools index ~{input_bam_10x_sorted}.subset.bam ~{input_bam_10x_sorted}.subset.bam.bai
    >>>

    output {
        File output_bam = "${input_bam_10x_sorted}.subset.bam"
        File output_bai = "${input_bam_10x_sorted}.subset.bam.bai"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Merge scRNA-seq and MAESTER bam files
task MergeBamFiles {
    input {
        File input_maester_bam
        File input_maester_bai
        File input_scrna_seq_bam
        File input_scrna_seq_bai
        String output_name = "merged_scRNAseq_and_MAESTER.bam"

        # Runtime
        String docker
        Int? disk_size_override
        Int cpu = 1
        Float memory = 3.5
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float input_scrna_seq_bam_size = size(input_scrna_seq_bam, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 3.0 * input_scrna_seq_bam_size)])

    command <<<
        set -euo pipefail

        samtools merge ~{output_name} ~{input_maester_bam} ~{input_scrna_seq_bam}
        samtools index ~{output_name} ~{output_name}.bai
    >>>

    output {
        File output_merged_bam = "${output_name}"
        File output_merged_bai = "${output_name}.bai"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Run maegatk with options -mr 3 to use UMIs with >=3 reads
# and -b to filter by high-quality CBs from scRNA-seq.
# Parallelized task - more cpu will speed up this process.
task Maegatk {
    input {
        File input_bai
        String output_name = "output.zip"

        # MAEGATK parameters
        File input_bam                   # -i; Input; a single, indexed bam file (required)
        String output_dir = "output_dir" # -o; Output directory for genotypes
        String? name                     # -n; Prefix for project name
        String? mito_genome              # -g; mitochondrial genome configuration. Requires bwa indexed fasta file or `rCRS` (built-in) [required]
        #String? ncores                  # -c; Number of cores to run the main job in parallel
        String? cluster                  # Message to send to Snakemake to execute jobs on cluster interface; see documentation.
        String? jobs                     # Max number of jobs to be running concurrently on the cluster interface
        String? barcode_tag              # -bt; Read tag (generally two letters) to separate single cells; valid and required only in 'bcall' mode
        File? barcodes                   # -b; File path to barcodes that will be extracted; useful only in 'bcall' mode.
        Int? min_barcode_reads           # -mb; Minimum number of mitochondrial reads for a barcode to be genotyped; useful only in 'bcall' mode; will not overwrite the --barcodes logic.
        Int? NHmax                       # Maximum number of read alignments allowed as governed by the NH flag. Default = 2
        Int? NMmax                       # Maximum number of paired mismatches allowed represented by the NM/nM tags. Default = 15
        Int? min_reads                   # -mr; Minimum number of supporting reads to call a consensus UMI/rread. Default = 1
        String? umi_barcode              # -ub; Read tag (generally two letters) to specify the UMI tag when removing duplicates for genotyping
        String? max_javamem              # -jm; Maximum memory for java for running duplicate removal. Default = 4000m
        Int? base_qual                   # -q; Minimum base quality for inclusion in the genotype count. Default = 0
        Int? alignment_quality           # -aq; Minimum alignment quality to include the read in genotype. Default = 0
        Int? nsamples                    # -ns; The number of samples/cells to be processed per iteration; default is all.
        String? keep_samples             # -k; Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file
        String? ignore_samples           # -x; Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file
        Boolean keep_temp_files = false  # -z; Keep all intermediate files
        Boolean skip_R = false           # -sr; Generate plain-text only output. Otherwise, this generates a .rds obejct that can be immediately read into R for downstream analysis.
        Boolean snake_stdout = true      # -so; Write snakemake log to stdout rather than a file.
        Boolean help = false             # Shows the help message and immediately exits.

        # MAEGATK mode
        Boolean bcall = false            # Required - bcall or support must be chosen to run maegatk.
        Boolean support = false

        # Runtime
        String docker
        Int? disk_size_override
        Int? cpu_override
        Float memory = 64.0
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float input_bam_size = size(input_bam, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 5.0 * input_bam_size)])
    Int cpu = select_first([cpu_override, if memory > 6.5 then 2 * floor(memory / 8) else 1])

	# --snake-stdout default to true: Attempts to fix "AttributeError in line 30 of /usr/local/lib/python3.6/dist-packages/maegatk/bin/snake/Snakefile.maegatk.Gather: 'InputFiles' object has no attribute 'depths'"
    # alias python=python3 - does not work - instead used a symbolic link in the Dockerfile. ln -s /usr/bin/python3 /usr/bin/python
    command <<<
        set -euo pipefail

        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8

        maegatk \
        --input=~{input_bam} \
        --ncores=~{cpu} \
        ~{if defined(output_dir) then "--output=~{output_dir} " else " "} \
        ~{if defined(name) then "--name=~{name} " else " "} \
        ~{if defined(mito_genome) then "--mito-genome=~{mito_genome} " else " "} \
        ~{if defined(cluster) then "--cluster=~{cluster} " else " "} \
        ~{if defined(jobs) then "--jobs=~{jobs} " else " "} \
        ~{if defined(barcode_tag) then "--barcode-tag=~{barcode_tag} " else " "} \
        ~{if defined(barcodes) then "--barcodes=~{barcodes} " else " "} \
        ~{if defined(min_barcode_reads) then "--min-barcode-reads=~{min_barcode_reads} " else " "} \
        ~{if defined(NHmax) then "--NHmax=~{NHmax} " else " "} \
        ~{if defined(NMmax) then "--NMmax=~{NMmax} " else " "} \
        ~{if defined(min_reads) then "--min-reads=~{min_reads} " else " "} \
        ~{if defined(umi_barcode) then "--umi-barcode=~{umi_barcode} " else " "} \
        ~{if defined(max_javamem) then "--max-javamem=~{max_javamem} " else " "} \
        ~{if defined(base_qual) then "--base-qual=~{base_qual} " else " "} \
        ~{if defined(alignment_quality) then "--alignment-quality=~{alignment_quality} " else " "} \
        ~{if defined(nsamples) then "--nsamples=~{nsamples} " else " "} \
        ~{if defined(keep_samples) then "--keep-samples=~{keep_samples} " else " "} \
        ~{if defined(ignore_samples) then "--ignore-samples=~{ignore_samples} " else " "} \
        ~{if keep_temp_files then "--keep-temp-files " else " "} \
        ~{if skip_R then "--skip-R " else " "} \
        ~{if snake_stdout then "--snake-stdout " else " "} \
        ~{if help then "--help " else " "} \
        ~{if bcall then "bcall " else " "} \
        ~{if support then "support " else " "}

        zip -r ~{output_name} ~{output_dir}
    >>>

    output {
        File output_zipped_file = "${output_name}"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

# Intersect common cell barcodes between scRNA-seq and MAESTER
# and assess coverage along mitochondrial genome.
task MtCoverage {
    input {
        File r_script           # 1.2_MT_Coverage - from the MAESTER github.
        File r_script_source    # 210215_FunctionsGeneral.R - from the MAESTER github.

        # Maester 1.2_MT_Coverage.R parameters
        File input_folder                                    # Zipped output from maegatk.
        String experiment_name                               # Sample name.
        String maegatk_full = "output_dir/final/maegatk.rds" # RDS file from maegatk.
        File metadata_df                                     # Metadata from scRNA-Seq.

        # Runtime
        String docker
        Int? disk_size_override
        Int cpu = 1
        Float memory = 20.0
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float input_folder_size = size(input_folder, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * input_folder_size)])
    String maegatk_full_current_dir = basename(maegatk_full)
    String metadata_df_current_dir = basename(metadata_df)

    command <<<
        set -euo pipefail
        unzip ~{input_folder}
        mv ~{maegatk_full} .
        mv ~{metadata_df} .

        R < ~{r_script} --no-save --args ~{r_script_source} ~{experiment_name} ~{maegatk_full_current_dir} ~{metadata_df_current_dir}
    >>>

    output {
        File output_plot = "${experiment_name}_plots.pdf"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}
