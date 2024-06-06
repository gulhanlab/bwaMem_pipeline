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

    call samtoolsSort {
        input:
            prefix = prefix,
            unsorted_bam = bwaMem.output_bam
    }

    call samtoolsIndex {
        input:
            prefix = prefix,
            sorted_bam = samtoolsSort.sorted_bam
    }

    call bamQC {
        input:
            prefix = prefix,
            bam_file = samtoolsSort.sorted_bam
    }

    output {
        #File aligned_bam = bwaMem.output_bam (remove; not needed, can just keep sorted_bam)
        File sorted_bam = samtoolsSort.sorted_bam
        File sorted_bam_idx = samtoolsIndex.sorted_bam_idx
        File genome_results = bamQC.genome_results
        File html_report = bamQC.html_report
        File pdf_report = bamQC.pdf_report
    }
}

# TASK DEFINITIONS #########################################
task bwaIndex {
    input {
        File reference_genome_fa
    }

    command {
        bwa index ${reference_genome_fa}
    }

    output {
        File bwa_index_amb = '${reference_genome_fa}.amb'
        File bwa_index_ann = '${reference_genome_fa}.ann'
        File bwa_index_bwt = '${reference_genome_fa}.bwt'
        File bwa_index_pac = '${reference_genome_fa}.pac'
        File bwa_index_sa = '${reference_genome_fa}.sa'
    }

    runtime {
        docker: 'gulhanlab/bwa-mem_pipeline:v1'
    }
}

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

    command {
        bwa mem ${reference_genome_fa} ${fastq1} ${fastq2} -t ${num_threads} > output.sam
        samtools view -S -b output.sam > ${prefix}_DNA.bam
    }
    
    output {
        File output_bam = '${prefix}_DNA.bam'
    }

    runtime {
        docker: 'gulhanlab/bwa-mem_pipeline:v1'
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}

task samtoolsSort {
    input {
        String prefix
        File unsorted_bam

        # Configurable
        Float memory
        Int disk_space
        Int num_threads
        Int num_preempt

        # For -m flag
        Int mem_per_thread = floor(memory / num_threads)
    }

    command {
        samtools sort -m ${mem_per_thread}G -@ ${num_threads} -o ${prefix}_DNA_sorted.bam ${unsorted_bam} 
    }

    output {
        File sorted_bam = '${prefix}_DNA_sorted.bam'
    }

    runtime {
        docker: 'gulhanlab/bwa-mem_pipeline:v1'
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}

task samtoolsIndex {
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
        samtools index -b -@ ${num_threads} -o ${prefix}_DNA_sorted.bam.bai ${sorted_bam}
    }

    output {
        File sorted_bam_idx = '${prefix}_DNA_sorted.bam.bai'
    }

    runtime {
        docker: 'gulhanlab/bwa-mem_pipeline:v1'
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}

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
        qualimap bamqc -bam ${bam_file} -outdir ${prefix}_DNA_qualimap_report -nt ${num_threads} --java-mem-size=${memory}G
    }

    output {
        File genome_results = '${prefix}_DNA_qualimap_report/genome_results.txt'
        File html_report = '${prefix}_qualimap_report/qualimapReport.html'
        File pdf_report = '${prefix}_qualimap_report/qualimapReport.pdf'
    }

    runtime {
        docker: 'gulhanlab/bwa-mem_pipeline:v1'
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
}
