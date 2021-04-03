nextflow run meripseqpipe \
-profile docker \
--genome GRCh37 \
--single_end 1 \
--readPaths 'SRP012098/SRR*/*fastq' \
--designfile 'metadata/meripmanifest.txt' \
--comparefile 'metadata/comparefile.txt' \
--aligners star \
--peakCalling_mode independence \
--peakMerged_mode rank \
--expression_analysis_mode DESeq2 \
--methylation_analysis_mode QNB \
--fasta refs/GRCh37.p13.genome.fa \
--gtf refs/gencode.v19.annotation.gtf \
--max_cpus 4 \
--max_memory '19GB'


--reads 'SRP012098/SRR*/*fastq' \