version 1.0

workflow bwaMem_pipeline {
    input {
        String prefix
        File fastq1
        File fastq2
        File reference_genome_fa
        # Optional
        File? reference_genome_amb
        File? reference_genome_ann
        File? reference_genome_bwt
        File? reference_genome_pac
        File? reference_genome_sa
    }

    # If BWA index files not provided, generate them
    if (!defined(reference_genome_amb) && !defined(reference_genome_ann) && !defined(reference_genome_bwt) &&
    !defined(reference_genome_pac) && !defined(reference_genome_sa)){
        call bwaIndex {
            input:
                reference_genome_fa=reference_genome_fa
        }
    }

    call bwaMem {
        input:
            prefix = prefix,
            fastq1 = fastq1,
            fastq2 = fastq2,
            reference_genome_fa = reference_genome_fa,
            reference_genome_amb = select_first([reference_genome_amb, bwaIndex.bwa_index_amb]),
            reference_genome_ann = select_first([reference_genome_ann, bwaIndex.bwa_index_ann]),
            reference_genome_bwt = select_first([reference_genome_bwt, bwaIndex.bwa_index_bwt]),
            reference_genome_pac = select_first([reference_genome_pac, bwaIndex.bwa_index_pac]),
            reference_genome_sa = select_first([reference_genome_sa, bwaIndex.bwa_index_sa])
    }

    call addOrReplaceReadGroups {
        input:
            prefix = prefix,
            sorted_bam = bwaMem.sorted_bam
    }

    # call samtoolsIndex {
    #     input:
    #         prefix = prefix,
    #         sorted_rg_bam = addOrReplaceReadGroups.sorted_rg_bam
    # }

    call bamQC {
        input:
            prefix = prefix,
            bam_file = addOrReplaceReadGroups.sorted_rg_bam
    }

    output {
        File sorted_rg_bam = addOrReplaceReadGroups.sorted_rg_bam
        File sorted_rg_bam_idx = samtoolsIndex.sorted_rg_bam_idx
        File genome_results = bamQC.genome_results
        File html_report = bamQC.html_report
    }
}

# TASK DEFINITIONS #########################################
task bwaIndex {
    input {
        File reference_genome_fa
    }

    command {
        bwa index ~{reference_genome_fa}
    }

    output {
        File bwa_index_amb = '~{reference_genome_fa}.amb'
        File bwa_index_ann = '~{reference_genome_fa}.ann'
        File bwa_index_bwt = '~{reference_genome_fa}.bwt'
        File bwa_index_pac = '~{reference_genome_fa}.pac'
        File bwa_index_sa = '~{reference_genome_fa}.sa'
    }

    runtime {
        docker: 'pegi3s/bwa:latest'
    }
}

# bwa mem and sort
task bwaMem {
    input {
        String prefix
        File fastq1
        File fastq2
        File reference_genome_fa
        File reference_genome_amb
        File reference_genome_ann
        File reference_genome_bwt
        File reference_genome_pac
        File reference_genome_sa

        # Configurable
        Float memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    # For -m flag for samtools sort
    Int mem_per_thread = floor(memory/(num_threads*2))

    command {
        bwa mem ~{reference_genome_fa} ~{fastq1} ~{fastq2} -t ~{num_threads} | \
        samtools sort -m ~{mem_per_thread}G -@ ~{num_threads} -o ~{prefix}_DNA_sorted.bam -
    }
    
    output {
        File sorted_bam = '~{prefix}_DNA_sorted.bam'
    }

    runtime {
        docker: 'gulhanlab/bwa_samtools:1.0'
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}

# Required for picard tools, generates bam index
task addOrReplaceReadGroups {
    input {
        String prefix
        File sorted_bam

        # Configurable
        Float memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command {
        java -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
            I=~{sorted_bam} \
            O=~{prefix}_DNA_sorted_rg.bam \
            RGID=1 \
            RGLB=lib_name \
            RGPL=ILLUMINA \
            RGPU=platform_unit \
            RGSM=~{prefix} \
            CREATE_INDEX=true
    }

    output {
        File sorted_rg_bam = '~{prefix}_DNA_sorted_rg.bam'
        File sorted_rg_bam_idx = '~{prefix}_DNA_sorted_rg.bam.bai'
    }

    runtime {
        docker: 'broadinstitute/picard:3.1.1'
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}

# task samtoolsIndex {
#     input {
#         String prefix
#         File sorted_rg_bam

#         # Configurable
#         Float memory
#         Int disk_space
#         Int num_threads
#         Int num_preempt
#     }

#     command {
#         samtools index -@ ~{num_threads} -o ~{prefix}_DNA_sorted_rg.bam.bai ~{sorted_rg_bam}
#     }

#     #samtools index -@ ~{num_threads} ~{sorted_rg_bam}

#     output {
#         #File sorted_rg_bam_idx = '~{sorted_rg_bam}.bai'
#         File sorted_rg_bam_idx = '~{prefix}_DNA_sorted_rg.bam.bai'
#     }

#     runtime {
#         docker: 'staphb/samtools:1.20'
#         memory: "~{memory} GB"
#         disks: "local-disk ~{disk_space} HDD"
#         cpu: "~{num_threads}"
#         preemptible: "~{num_preempt}"
#     }
# }

task bamQC {
    input {
        String prefix
        File bam_file

        # Configurable
        Int memory # java-mem-size must be int
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command {
        qualimap bamqc \
            -bam ~{bam_file} \
            -outdir ~{prefix}_DNA_qualimap_report \
            -nt ~{num_threads} \
            --java-mem-size=~{memory}G
    }

    output {
        File genome_results = '~{prefix}_DNA_qualimap_report/genome_results.txt'
        File html_report = '~{prefix}_DNA_qualimap_report/qualimapReport.html'
    }

    runtime {
        docker: 'pegi3s/qualimap:2.2.1'
        memory: "~{memory} GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: "~{num_preempt}"
    }
}
