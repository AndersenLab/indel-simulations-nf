

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
    	set file(indel_positions), file(simulated_genome), indel, size, file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') into sim_index

    
    """
    	bwa index -a bwtsw ${simulated_genome}
    """
}


process perform_alignment {

	echo true

    input:
        set file(indel_positions), file(simulated_genome), indel, size,file('simulated_genome.fa.amb'), file('simulated_genome.fa.ann'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') from sim_index
    output:
    	set file(indel_positions), file(simulated_genome), indel, size, file('ECA316.bam'), file('ECA316.bam.bai'), file('simulated_genome.fa.bwt'), file('simulated_genome.fa.pac'), file('simulated_genome.fa.sa') into fq_to_bam

    
    """
    	mkdir temp
        bwa mem -t ${alignment_cores} -R '@RG\tID:N2\tLB:N2\tSM:N2' ${simulated_genome} ${input_fastq1} ${input_fastq2} | \\
        sambamba view --nthreads=${alignment_cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${alignment_cores} --show-progress --tmpdir=temp --out=${sample_id}.bam /dev/stdin
        sambamba index --nthreads=${alignment_cores} ${sample_id}.bam
    """
}







