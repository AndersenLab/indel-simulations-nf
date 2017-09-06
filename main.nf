

//input_fasta = file("/projects/b1059/data/genomes/c_elegans/WS245/WS245.tmp.fa")
input_fasta = "/Users/stefan/Indels/c_elegans.PRJNA13758.WS245.genomic.fa"
input_fastq1 = "/Users/stefan/AndersenLab/Github_Repos/indel-simulations-nf/data/ECA316_1P.fq.gz"
input_fastq2 = "/Users/stefan/AndersenLab/Github_Repos/indel-simulations-nf/data/ECA316_2P.fq.gz"
alignment_cores = 4
sample_id = "ECA316"

//indel = Channel.from(['insertion', 'deletion'])
//size = Channel.from(['1,10', '11,100','101,1000','1001,20000'])
indel = Channel.from(['deletion'])
size = Channel.from(['1,10'])

indel.spread(size).into { out_use; out_print}

//out_print.subscribe { println it }

process rsvsim {

	echo true

	tag { "${size} - ${indel}" }

    input:
       set indel, size from out_use
  
    output:
    	set file('indel_positions.tsv'), file('simulated_genome.fa'), indel, size into variant_positions

    """
    #!/usr/bin/env Rscript --vanilla

    library(readr)
    library(RSVSim)
	library(Biostrings)

    range <- strsplit("${size}", ",")

	small_bound <- as.numeric(range[[1]][1])
	large_bound <- as.numeric(range[[1]][2])

	if("${indel}" == "insertion"){
		simulated_genome <- simulateSV(output = NA, genome = "${input_fasta}" , 
                   dels = 10, 
                   sizeDels = sample(seq(small_bound,large_bound, by = 1), size = 1000, replace = T))

                   readr::write_tsv(metadata(simulated_genome)[[1]], 'indel_positions.tsv')
                   Biostrings::writeXStringSet(simulated_genome,'simulated_genome.fa')

	} else {
		simulated_genome <- simulateSV(output = NA, genome = "${input_fasta}", 
                   dels = 10, 
                   sizeDels = sample(seq(small_bound,large_bound, by = 1), size = 1000, replace = T)) 

                   readr::write_tsv(metadata(simulated_genome)[[1]], 'indel_positions.tsv')
                   Biostrings::writeXStringSet(simulated_genome,'simulated_genome.fa')
	}

    """
}

process index_simulated_genome {

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size from variant_positions
    output:
    	set file(indel_positions), file(simulated_genome), indel, size, file('simulated_genome.fa.fai'), file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') into sim_index

    
    """
    	bwa index -a bwtsw ${simulated_genome}
    	samtools faidx ${simulated_genome}
    """
}


process perform_alignment {

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('simulated_genome.fa.fai'),file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from sim_index
    output:
    	set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') into fq_to_bam

    
    """
    	mkdir temp
        bwa mem -t ${alignment_cores} -R '@RG\tID:N2\tLB:N2\tSM:N2' ${simulated_genome} ${input_fastq1} ${input_fastq2} | \\
        sambamba view --nthreads=${alignment_cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${alignment_cores} --show-progress --tmpdir=temp --out=${sample_id}.bam /dev/stdin
        sambamba index --nthreads=${alignment_cores} ${sample_id}.bam
    """
}

//split fq_to_bam into multiple channels for variant calling
fq_to_bam.into { platypus_files; freebayes_files; gatk_files; bcftools_files; pindel_files; delly_files; lumpy_files; post_caller_process_files }

process platypus {

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from platypus_files
    output:
    	set file('platypus_output.vcf') into platypus_vcf


	    """
	    python ~/Indels/Platypus_0.8.1/Platypus.py callVariants --bamFiles=ECA316.bam --refFile=${simulated_genome} --output=platypus_output.vcf

    	"""
}


process freebayes {

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from freebayes_files
    output:
    	set file('freebayes_output.vcf') into freebayes_vcf


	    """
	    freebayes --fasta-reference simulated_genome.fa ECA316.bam | vcffilter -f "QUAL > 20" > freebayes_output.vcf 
    	"""
}

process gatk {

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from gatk_files
    output:
    	set file('gatk_output.vcf') into gatk_vcf


	    """
	    picard CreateSequenceDictionary R=simulated_genome.fa O=simulated_genome.dict
	    java -jar /Users/stefan/Indels/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
     	-R simulated_genome.fa \
     	-T HaplotypeCaller \
     	-I ECA316.bam \
     	--emitRefConfidence GVCF \
     	-o gatk_output.vcf \
     	-variant_index_type LINEAR \
     	-variant_index_parameter 128000
    	"""
}

process bcftools {
		echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from bcftools_files
    output:
    	set file('bcftools_output.vcf') into bcftools_vcf


	    """
		bcftools mpileup -f simulated_genome.fa ECA316.bam | bcftools call -mv -Ov -o bcftools_output.bcf
		bcftools convert -Ov -o bcftools_output.vcf bcftools_output.bcf
    	"""
}

process delly {
		echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from delly_files
    output:
    	set file('delly_output.vcf') into delly_vcf

    	script:
    	if( indel == 'deletion')
	    """
		delly call -t DEL -g simulated_genome.fa -o delly_output.bcf ECA316.bam
		bcftools convert -Ov -o delly_output.vcf delly_output.bcf
    	"""
    	else
    	"""
		delly call -t DEL -g simulated_genome.fa -o delly_output.bcf ECA316.bam
		bcftools convert -Ov -o delly_output.vcf delly_output.bcf
    	"""

}

process lumpy{

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from lumpy_files
    output:
    	set file('lumpy_output.vcf') into lumpy_vcf
	    

	    """
		samtools view -b -F 1294 ECA316.bam > ECA316.discordants.unsorted.bam
		samtools view -h ECA316.bam \
    	| /Users/stefan/AndersenLab/Github_Repos/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    	| samtools view -Sb - \
    	> ECA316.splitters.unsorted.bam

    	samtools sort -o ECA316.discordants.bam ECA316.discordants.unsorted.bam 
		samtools sort -o ECA316.splitters.bam ECA316.splitters.unsorted.bam 

		lumpyexpress \
    	-B ECA316.bam \
    	-S ECA316.splitters.bam \
    	-D ECA316.discordants.bam \
    	-o lumpy_output.vcf
    	"""

}

//process pindel {
//	echo true
//
//   input:
//        set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.fai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from pindel_files
//    output:
//    	set file('pindel_output.vcf') into pindel_vcf
//
//
//    	"""
//    	echo 'ECA316.bam     300   pindel_output' >pindelconfig.txt
//    	link="\$(readlink simulated_genome.fa)"
//    	rm simulated_genome.fa
//    	ls -la
//    	cp "\${link}" .
//    	docker run opengenomics/pindel pindel -f ./simulated_genome.fa -i ./pindelconfig.txt -c ALL -o ./pindel_output
//    	docker run opengenomics/pindel pindel2vcf -P ./pindel_output -r ./simulated_genome.fa -R simulated -d 2017
//    	"""
//}

process generate_tsv {

	echo true

    input:
    	set file('platypus_output.vcf') from platypus_vcf
    	set file('freebayes_output.vcf') from freebayes_vcf
    	set file('gatk_output.vcf') from gatk_vcf
    	set file('bcftools_output.vcf') from bcftools_vcf
    	set file('delly_output.vcf') from delly_vcf
        set file('lumpy_output.vcf') from lumpy_vcf

    output:
    	set file('platypus_output.tsv'), file('freebayes_output.tsv'), file('GATK_output.tsv'), file('bcftools_output.tsv'), file('delly_output.tsv'), file('lumpy_output.tsv') into all_tsv_outputs


	    """
	    vcf2tsv platypus_output.vcf -g > platypus_output.tsv
	    vcf2tsv freebayes_output.vcf -g > freebayes_output.tsv
	    vcf2tsv GATK_output.vcf -g > GATK_output.tsv
	    vcf2tsv bcftools_output.vcf -g > bcftools_output.tsv
		vcf2tsv delly_output.vcf -g > delly_output.tsv
		vcf2tsv lumpy_output.vcf -g > lumpy_output.tsv
    	"""


}
