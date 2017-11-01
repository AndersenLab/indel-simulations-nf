#!/usr/bin/env nextflow 
fastqdir = config.fastqdir
input_fasta = config.input_fasta
alignment_cores = config.alignment_cores
sample_id = config.sample_id
sambamba = config.sambamba
platypus = config.platypus
lumpy_extract_split = config.lumpy_extract_split
gatkdir = config.gatkdir
picarddir = config.picarddir
indel_simulator = config.indel_simulator

// Define ranges of indel sizes
indel = Channel.from(['insertion', 'deletion'])
size = Channel.from(['1 10', '11 100','101 1000','1001 20000'])

indel.spread(size).into { out_use; out_print}

// Initialize channels to input fastq files
id = Channel.create()
fq1 = Channel.create()
fq2 = Channel.create()

// load FASTQ files into fastq_pairs
Channel
        .fromFilePairs(fastqdir, flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: ${fastqdir}" }
        .separate( id, fq1, fq2)


// combine fastq files 
process conc {
     
    echo true

    input:
    val fqs1 from fq1.toSortedList()
    val fqs2 from fq2.toSortedList()

    output:
    set file('fq1.fq.gz'), file('fq2.fq.gz') into fq_pairs_combined
 
    script:
    """
    zcat ${fqs1.join(' ')} | gzip > fq1.fq.gz
    zcat ${fqs2.join(' ')} | gzip > fq2.fq.gz
    """
}

// simulate indels
process simulate {

    //publishDir "results/${indel}/${size}", mode: 'copy'

    tag { "${size} - ${indel}" }

    input:
        set indel, size from out_use
  
    output:
        set file('indel_positions.tsv'), file('simulated_genome.fa'), indel, size into variant_positions

    "python ${indel_simulator} -g ${input_fasta} -i $indel -r $size -n 1000 -l 200"
}

// index simulated genome
process index_simulated_genome {

    echo true

    input:
        set file(indel_positions), file('simulated_genome.fa'), indel, size from variant_positions
    output:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') into indexed_simulated_genome

    
    """
        bwa index -a bwtsw simulated_genome.fa
        samtools faidx simulated_genome.fa
        tail -n +2 indel_positions.tsv | sort -k 1,1 -k 2,2n | grep -vwE MtDNA > indel_positions.bed
    """
}

// align reads to simulated genome
process perform_alignment {

    echo true

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from indexed_simulated_genome
        set file('fq1.fq.gz'), file('fq2.fq.gz') from fq_pairs_combined
    output:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai') into aligned_bam_bai
    
    """
        bwa mem -t ${alignment_cores} -R '@RG\tID:N2\tLB:N2\tSM:N2\tPL:ILLUMINA' simulated_genome.fa fq1.fq.gz fq2.fq.gz | \\
        ${sambamba} view --nthreads=${alignment_cores} --sam-input --format=bam --with-header /dev/stdin | \\
        ${sambamba} sort --nthreads=${alignment_cores} --show-progress --tmpdir=temp --out=${sample_id}.bam /dev/stdin
        ${sambamba} index --nthreads=${alignment_cores} ${sample_id}.bam
    """
}




process gatk_base_recalibrate {

    echo true

    publishDir "results/${indel}/${size}/plots", mode: 'copy', pattern: '*.pdf'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai') from aligned_bam_bai
    output:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') into gatk_recalibrated
        file('recalibration_plots.pdf') into gatk_recalibrated_images

        """
        java -jar ${picarddir} CreateSequenceDictionary R=simulated_genome.fa O=simulated_genome.dict

        java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \
            -T BaseRecalibrator \
            -R simulated_genome.fa \
            -I N2.bam \
            -L I:2000000-12000000 \
            -knownSites indel_positions.bed \
            -o recal_data.table 

        java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \
            -T BaseRecalibrator \
            -R simulated_genome.fa \
            -I N2.bam \
            -L I:2000000-12000000 \
            -knownSites indel_positions.bed \
            -BQSR recal_data.table  \
            -o post_recal_data.table 

        java -jar ${gatkdir} \
            -T AnalyzeCovariates \
            -R simulated_genome.fa \
            -L I:2000000-12000000 \
            -before recal_data.table \
            -after post_recal_data.table \
            -plots recalibration_plots.pdf

        java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \
            -T PrintReads \
            -R simulated_genome.fa \
            -I N2.bam \
            -BQSR recal_data.table \
            --defaultBaseQualities 0 \
            -o N2_recald.bam 
        """
}

gatk_recalibrated.into {gatk_recalibrated_gatk;
                        gatk_recalibrated_platypus;
                        gatk_recalibrated_freebayes;
                        gatk_recalibrated_gridss;
                        gatk_recalibrated_delly;
                        gatk_recalibrated_lumpy;
                        gatk_recalibrated_bcftools;
                        gatk_recalibrated_post_caller1;
                        gatk_recalibrated_post_caller2}



process gatk_call {

    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_gatk
    output:
        set file("gatk_output_indel.vcf"), file("gatk_output_indel.vcf.idx"), file("gatk_window.tsv"), file("gatk_output_indel.tsv") into gatk_vcf

        """
        java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \
            -T HaplotypeCaller \
            -R simulated_genome.fa \
            -I N2_recald.bam \
            --genotyping_mode DISCOVERY \
            -stand_call_conf 30 \
            --sample_ploidy 1 \
            -o gatk_output.vcf

        java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \
           -R simulated_genome.fa \
           -T SelectVariants \
           -V gatk_output.vcf \
           -o gatk_output_indel.vcf \
           -selectType INDEL



        convert2bed --input=vcf --output=bed <gatk_output_indel.vcf> gatk_output.bed

        bedtools window -a indel_positions.bed -b gatk_output.bed -w 100 > gatk_window.tsv

        vk vcf2tsv wide --print-header gatk_output_indel.vcf > gatk_output_indel.tsv

        """
}

process platypus_call {

    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_platypus
    output:
        set file('platypus_output.vcf'), file("platypus_window.tsv"), file("platypus_output.tsv") into platypus_vcf


        """
        python ${platypus} callVariants --nCPU=${alignment_cores} --bamFiles=N2_recald.bam --assemble=1 --refFile=simulated_genome.fa --output=platypus_output.vcf

        convert2bed --input=vcf --output=bed <platypus_output.vcf> platypus_output.bed

        bedtools window -a indel_positions.bed -b platypus_output.bed -w 100 > platypus_window.tsv

        vk vcf2tsv wide --print-header platypus_output.vcf > platypus_output.tsv
        """
}

process freebayes_call {

    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_freebayes
    output:
        set file('freebayes_output.vcf'), file("freebayes_window.tsv"), file("freebayes_output.tsv") into freebayes_vcf


        """
        freebayes --fasta-reference simulated_genome.fa --ploidy 1 N2_recald.bam > freebayes_output.vcf 

        convert2bed --input=vcf --output=bed <freebayes_output.vcf> freebayes_output.bed

        bedtools window -a indel_positions.bed -b freebayes_output.bed -w 100 > freebayes_window.tsv

        vk vcf2tsv wide --print-header freebayes_output.vcf > freebayes_output.tsv
        """
}

process gridss_call {

    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_gridss
    output:
        set file('gridss_output.sv.vcf'), file("gridss_window.tsv"), file("gridss_output.tsv") into gridss_vcf
        

        """
        java -ea -Xmx31g -Dsamjdk.create_index=true \
            -Dsamjdk.use_async_io_read_samtools=true \
            -Dsamjdk.use_async_io_write_samtools=true \
            -Dsamjdk.use_async_io_write_tribble=true \
            -cp /projects/b1059/software/gridss-1.4.3-jar-with-dependencies.jar gridss.CallVariants \
            REFERENCE_SEQUENCE=simulated_genome.fa \
            TMP_DIR=. \
            WORKING_DIR=. \
            INPUT=N2_recald.bam \
            OUTPUT=gridss_output.sv.vcf \
            ASSEMBLY=gridss_output.gridss.assembly.bam

        convert2bed --input=vcf --output=bed <gridss_output.sv.vcf> gridss_output.bed

        bedtools window -a indel_positions.bed -b gridss_output.bed -w 100 > gridss_window.tsv

        vk vcf2tsv wide --print-header gridss_output.sv.vcf > gridss_output.tsv        
        """

}

process delly_call {
        
    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_delly
    output:
        set file('delly_output.vcf'), file("delly_window.tsv"), file("delly_output.tsv") into delly_vcf

        script:
        if( indel == 'deletion')
        """
        delly call -t INS -g simulated_genome.fa -o delly_output.bcf N2_recald.bam
        bcftools convert -Ov -o delly_output.vcf delly_output.bcf

        convert2bed --input=vcf --output=bed <delly_output.vcf> delly_output.bed

        bedtools window -a indel_positions.bed -b delly_output.bed -w 100 > delly_window.tsv

        vk vcf2tsv wide --print-header delly_output.vcf > delly_output.tsv 
        """
        else
        """
        delly call -t DEL -g simulated_genome.fa -o delly_output.bcf N2_recald.bam
        bcftools convert -Ov -o delly_output.vcf delly_output.bcf

        convert2bed --input=vcf --output=bed <delly_output.vcf> delly_output.bed

        bedtools window -a indel_positions.bed -b delly_output.bed -w 100 > delly_window.tsv

        vk vcf2tsv wide --print-header delly_output.vcf > delly_output.tsv 
        """

}

process lumpy_call {

    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_lumpy
    output:
        set file('lumpy_output.vcf'), file("lumpy_window.tsv"), file("lumpy_output.tsv") into lumpy_vcf
        

        """
        samtools view -b -F 1294 N2_recald.bam > N2.discordants.unsorted.bam
        samtools view -h N2_recald.bam \
        | ${lumpy_extract_split} -i stdin \
        | samtools view -Sb - \
        > N2.splitters.unsorted.bam

        samtools sort N2.discordants.unsorted.bam -O bam -T N2.discordants -o N2.discordants.bam  
        samtools sort N2.splitters.unsorted.bam -O bam -T N2.splitters -o N2.splitters.bam 

        lumpyexpress \
        -B N2_recald.bam \
        -S N2.splitters.bam \
        -D N2.discordants.bam \
        -o lumpy_output.vcf

        convert2bed --input=vcf --output=bed <lumpy_output.vcf> lumpy_output.bed

        bedtools window -a indel_positions.bed -b lumpy_output.bed -w 100 > lumpy_window.tsv

        vk vcf2tsv wide --print-header lumpy_output.vcf > lumpy_output.tsv

        """

}

process bcftools_call {
    
    echo true

    publishDir "results/${indel}/${size}", mode: 'copy'

    input:
        set file('indel_positions.bed'), indel, size, file("simulated_genome.fa"), file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa'), file('N2.bam'), file('N2.bam.bai'), file('N2_recald.bam'), file('N2_recald.bai'), file('simulated_genome.dict') from gatk_recalibrated_bcftools
    output:
        set file('bcftools_output.vcf'), file("bcftools_window.tsv"), file("bcftools_output.tsv") into bcftools_vcf


        """
        echo -e N2 1 > sample_file.ped

        bcftools mpileup --threads ${alignment_cores} -f simulated_genome.fa N2_recald.bam | bcftools call -mv -Ov -o bcftools_output.bcf
        bcftools convert -Ov -o bcftools_output.vcf bcftools_output.bcf

        convert2bed --input=vcf --output=bed <bcftools_output.vcf> bcftools_output.bed

        bedtools window -a indel_positions.bed -b bcftools_output.bed -w 100 > bcftools_window.tsv

        vk vcf2tsv wide --print-header bcftools_output.vcf > bcftools_output.tsv
        """
}

