from glob import glob

#mappers = ['fairy', 'minimap', 'bwa']

reads =  "QC/reads/"
#reads =  "reads/"

PAIR1 = "_R1.fastq.gz"
PAIR2 = "_R2.fastq.gz"

#PAIR1 = "_1.fastq.gz"
#PAIR2 = "_2.fastq.gz"

SAMPLETYPE = ['all']
BINNERS = ['maxbin']
#BINNERS = ['metabat2']
####Metabinner
CURR_DIR = "chicken/"
cc_covfiles = "concoct_covfiles/"

#TAG = "SRR"
#TAG = "ERR"
TAG = "ERR"

ASS = "/assembly/scaffolds.fasta"
#ASS = "/contigs.fasta.gz"
#ASS = "/assembly.fasta"

covfiles = "covfiles/"
bins = "bins/"
benchmarks = "benchmarks/"
envs = "/mnt/disks/MOUNT_DIR/envs/"
checkm2 = "checkm2_results/"

#reads1=reads + TAG + "{sample}" + PAIR1,
#reads2=reads + TAG + "{sample}" + PAIR2,

SAMPLES = glob_wildcards(reads + TAG + "{sample}" + PAIR1).sample
#SAMPLES = glob_wildcards(reads + TAG + "{sample}.fastq.gz").sample
#SAMPLES = glob_wildcards(TAG + "{sample}").sample
#SAMPLES = ["19064410", "19064411"]

long_threads = 40
short_threads = 20

SEQ_TYPE = ['short']
if SEQ_TYPE[0] == 'long':
    MAPPERS = ['fairy', 'minimap']
    long_reads = expand(reads + TAG + "{sample}.fastq.gz", sample= SAMPLES)
    reads1 = []
    reads2 = []
    long_sketches = expand("fairy_sketches_long/" + TAG + "{sample}" + "fastq.gz.bcsp", sample=SAMPLES)
    short_sketches = []
else:
    MAPPERS = ['fairy', 'bwa']
    long_reads = []
    reads1=expand(reads + TAG + "{sample}" + PAIR1, sample = SAMPLES)
    print(reads1)
    reads2=expand(reads + TAG + "{sample}" + PAIR2, sample = SAMPLES)
    print(reads2)
    short_sketches = expand("fairy_sketches_short/" + TAG + "{sample}" + PAIR1 + ".paired.bcsp", sample=SAMPLES)
    long_sketches = []

rule all:
    input:
        expand(checkm2 + "{TAG}{sample}_{mapper}_{sampletype}_{binner}_{seqtype}_cm2", TAG=TAG, sample=SAMPLES, mapper=MAPPERS,sampletype=SAMPLETYPE,binner=BINNERS,seqtype=SEQ_TYPE)
        #expand(checkm2 + "{TAG}{sample}_{mapper}_{sampletype}_{binner}_long_cm2", TAG=TAG, sample=SAMPLES, mapper=MAPPERS,sampletype=SAMPLETYPE,binner=BINNERS)
        #expand("{TAG}{sample}_minimap_single_covfile.tsv", TAG=TAG, sample=SAMPLES)


rule coverm_minimap_pair_single:
    input:
        reads1=reads + TAG + "{sample}" + PAIR1,
        reads2=reads + TAG + "{sample}" + PAIR2,
        fasta=TAG + "{sample}" + ASS
    output:
        cov_file=covfiles +  TAG + "{sample}_minimap_single_short.tsv"
    threads:
        long_threads
    benchmark:
        benchmarks + TAG + "{sample}_minimap_all_benchmark.txt"
    shell:
        "coverm contig -t {threads} -m metabat -1 {input.reads1} -2 {input.reads2} --reference {input.fasta} > {output}"


rule coverm_minimap_pair_all:
    input:
        reads1=reads1,
        reads2=reads2,
        fasta=TAG + "{sample}" + ASS
    output:
        cov_file=covfiles +  TAG + "{sample}_minimap_all_short.tsv"
    threads:
        long_threads
    benchmark:
        benchmarks + TAG + "{sample}_minimap_all_benchmark.txt"
    shell:
        "coverm contig -t {threads} -m metabat -1 {input.reads1} -2 {input.reads2} --reference {input.fasta} > {output}"


rule coverage_minimap2_single:
    input:
        fasta=TAG + "{sample}" + ASS,
        read=reads +TAG + "{sample}.fastq.gz" 
    output:
        cov_file=covfiles + TAG + "{sample}_minimap_single_long.tsv"
    threads:
        long_threads
    benchmark:
        benchmarks + TAG + "{sample}_minimap_single_benchmark.txt"
    run:
        shell("minimap2 -a {input.fasta} {input.read} -t {threads} | samtools sort -@ {threads} -o {wildcards.sample}.bam")
        shell("jgi_summarize_bam_contig_depths --outputDepth {output.cov_file} *.bam")
        shell("rm *.bam")


rule coverage_minimap2:
   input:
       fasta=TAG + "{sample}" + ASS,
       read=long_reads
   output:
       cov_file=covfiles + TAG + "{sample}_minimap_all_long.tsv"
   threads:
       long_threads
   benchmark:
       "benchmarks/" + TAG + "{sample}_minimap_all_benchmark.txt"
   run:
       for i,f in enumerate(input.read):
           shell("minimap2 -a {input.fasta} {f} -t {threads} | samtools sort -@ {threads} -o {i}.bam")
       shell("jgi_summarize_bam_contig_depths --outputDepth {output.cov_file} *.bam")
       shell("rm *.bam")

rule coverm_strobealign_single:
    input:
        reads1=reads + TAG + "{sample}" + PAIR1,
        reads2=reads + TAG + "{sample}" + PAIR2,
        fasta=TAG + "{sample}" + ASS
    output:
        cov_file=covfiles +  TAG + "{sample}_strobealign_single_short.tsv"
    threads:
        short_threads
    benchmark:
        benchmarks + TAG + "{sample}_strobealign_single_benchmark.txt"
    shell:
        "coverm contig -p strobealign -t {threads} -m metabat -1 {input.reads1} -2 {input.reads2} --reference {input.fasta} > {output}"


rule coverm_strobealign_all:
    input:
        reads1=reads1,
        reads2=reads2,
        fasta=TAG + "{sample}" + ASS
    output:
        cov_file=covfiles +  TAG + "{sample}_strobealign_all_short.tsv"
    threads:
        short_threads
    benchmark:
        benchmarks + TAG + "{sample}_strobealign_all_benchmark.txt"
    shell:
        "coverm contig -p strobealign -t {threads} -m metabat -1 {input.reads1} -2 {input.reads2} --reference {input.fasta} > {output}"


rule coverm_bwa_single:
    input:
        reads1=reads + TAG + "{sample}" + PAIR1,
        reads2=reads + TAG + "{sample}" + PAIR2,
        fasta=TAG + "{sample}" + ASS
    output:
        cov_file=covfiles +  TAG + "{sample}_bwa_single_short.tsv"
    threads:
        short_threads
    benchmark:
        benchmarks + TAG + "{sample}_bwa_single_benchmark.txt"
    shell:
        "coverm contig -p bwa-mem -t {threads} -m metabat -1 {input.reads1} -2 {input.reads2} --reference {input.fasta} > {output}"


rule coverm_bwa_all:
    input:
        reads1=reads1,
        reads2=reads2,
        fasta=TAG + "{sample}" + ASS
    output:
        cov_file=covfiles +  TAG + "{sample}_bwa_all_short.tsv"
    threads:
        short_threads
    benchmark:
        benchmarks + TAG + "{sample}_bwa_all_benchmark.txt"
    shell:
        "coverm contig -p bwa-mem -t {threads} -m metabat -1 {input.reads1} -2 {input.reads2} --reference {input.fasta} > {output}"

rule fairy_sketch_short:
    input:
        reads1=reads1,
        reads2=reads2
    output:
        sketches = short_sketches
    params:
        sdir = "fairy_sketches_short"
    threads:
        short_threads
    benchmark:
        "benchmarks/fairy_sketch_time.txt"
    shell:
        "fairy sketch --debug -1 {input.reads1} -2 {input.reads2} -d {params.sdir} -t {threads}"

rule fairy_sketch_long:
    input:
        long_reads = long_reads
    output:
        sketches = long_sketches
    params:
        sdir = "fairy_sketches_long"
    threads:
        short_threads
    benchmark:
        "benchmarks/fairy_sketch_time.txt"
    shell:
        "fairy sketch -r {input.long_reads} -d {params.sdir} -t {threads}"

rule fairy_profile_1_long:
    input:
        fasta=TAG + "{sample}" + ASS,
        sketches = "fairy_sketches_long/" + TAG + "{sample}" + "fastq.gz.bcsp"
    output:
        tsv=covfiles + TAG + "{sample}_fairy_single_long.tsv"
    threads:
        short_threads
    shell:
        "fairy coverage {input.sketches} {input.fasta} -t {threads} -o {output.tsv}"


rule fairy_profile_long:
    input:
        fasta=TAG + "{sample}" + ASS,
        sketches = long_sketches
    output:
        tsv=covfiles + TAG + "{sample}_fairy_all_long.tsv"
    threads: 
        short_threads
    shell:
        "fairy coverage {input.sketches} {input.fasta} -t {threads} -o {output.tsv}"

rule fairy_profile_1_short:
    input:
        fasta=TAG + "{sample}" + ASS,
        sketches = "fairy_sketches_short/" + TAG + "{sample}" + PAIR1 + ".paired.bcsp"
    output:
        tsv=covfiles + TAG + "{sample}_fairy_single_short.tsv"
    threads:
        short_threads
    shell:
        "fairy coverage {input.sketches} {input.fasta} -t {threads} -o {output.tsv}"

rule fairy_profile_short:
    input:
        fasta=TAG + "{sample}" + ASS,
        sketches = short_sketches
    output:
        tsv=covfiles + TAG + "{sample}_fairy_all_short.tsv"
    threads: 
        short_threads
    shell:
        "fairy coverage {input.sketches} {input.fasta} -t {threads} -o {output.tsv}"

#rule metabat2_fairy1:
#    input:
#        fasta="ERR{sample}/contigs.fasta.gz",
#        sorted="ERR{sample}-fairy1.tsv"
#    output:
#        bin=directory("ERR{sample}_fairy_Single_mags")
#    params:
#        out="ERR{sample}_fairy_Single_bin"
#    threads:
#        25
#    shell:
#        """
#        metabat2 -i {input.fasta} -a {input.sorted} -t 25 --outFile {params.out} --seed 13 --minContig 1500
#        mkdir {output.bin}
#        mv {params.out}*.fa {output.bin}
#        """
#
#
rule metabat2:
    input:
        fasta=TAG + "{sample}" + ASS,
        sorted= covfiles + TAG + "{sample}_{mapper}_{sampletype}_{seqtype}.tsv"
    output:
        bin= directory(bins + TAG + "{sample}_{mapper}_{sampletype}_metabat2_{seqtype}")
    params:
        out=TAG + "{sample}_{mapper}_{sampletype}_metabat2_{seqtype}_bin"
    threads:
        25
    shell:
        """
        metabat2 -i {input.fasta} -a {input.sorted} -t {threads} --outFile {params.out} --seed 13 --minContig 1500
        mkdir -p {output.bin}
        mv {params.out}*.fa {output.bin}
        """

rule vamb:
    input:
        fasta=TAG + "{sample}" + ASS,
        covfile = covfiles + TAG + "{sample}_{mapper}_{sampletype}_{seqtype}.tsv"
    output:
        out= directory(bins + TAG + "{sample}_{mapper}_{sampletype}_vamb_{seqtype}")
    threads:
        short_threads
    conda:
        envs + "vamb_env.yaml"
    shell:
        """
        vamb --outdir {output.out} --jgi {input.covfile} --fasta {input.fasta} --minfasta 200000
        mv {output.out}/bins/* {output.out}/
        """

rule checkm2_predict_mb2:
    input:
        folder = bins + "{TAG}{sample}_{mapper}_{sampletype}_metabat2_{seqtype}"
    output:
        cm2 = directory(checkm2 + "{TAG}{sample}_{mapper}_{sampletype}_metabat2_{seqtype}_cm2")
    conda:
        envs + "checkm2_env.yaml"
    threads: short_threads # Adjust the number of threads as necessary
    shell:
        "~/checkm2/bin/checkm2 predict --threads {threads} --input {input.folder} --output-dir {output.cm2} -x fa --force;"
        "rm {output.cm2}/protein_files/*;"
        "rm {output.cm2}/diamond*/*"

rule checkm2_predict_vamb:
    input:
        folder = bins + "{TAG}{sample}_{mapper}_{sampletype}_vamb_{seqtype}"
    output:
        cm2 = directory(checkm2 + "{TAG}{sample}_{mapper}_{sampletype}_vamb_{seqtype}_cm2")
    conda:
        envs + "checkm2_env.yaml"
    threads: long_threads
    shell:
        "~/checkm2/bin/checkm2 predict --threads {threads} --input {input.folder}/ --output-dir {output.cm2} -x fna --force;"
        "rm {output.cm2}/protein_files/*;"
        "rm {output.cm2}/diamond*/*"

rule trim_mapping_cc:
    input: 
        sorted= covfiles + TAG + "{sample}_"+"{mapper}_all_short.tsv",
    output:
        trim = cc_covfiles +  TAG + "{sample}_"+"{mapper}_all_short_trim.tsv",
    run: 
        import csv

        # Define the input and output files
        input_file = input['sorted']
        output_file = output['trim']

        def filter_rows(input_file, output_file):
            with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
                # Create reader and writer objects
                reader = csv.reader(infile, delimiter='\t')
                writer = csv.writer(outfile, delimiter='\t')

                # Read and write the header
                header = next(reader)
                print(header)
                abridged_header = [header[0]] + [header[i] for i in range(3,len(header),2)]
                writer.writerow(abridged_header)

                # Iterate over the rows
                for row in reader:
                    abridged_row = [row[0]] + [row[i] for i in range(3,len(row),2)]
                    # Split the 'Contig_name' field and extract the length
                    length = int(row[1])
                    # Check if the length is greater than 1500
                    if length > 1500:
                        writer.writerow(abridged_row)

        # Call the function with the file paths
        print(input_file, output_file)
        filter_rows(input_file, output_file)


rule metabinner_fairy:
    input:
        fasta=TAG + "{sample}" + ASS,
        trim = cc_covfiles +  TAG + "{sample}_"+"{mapper}_all_short_trim.tsv"
    output:
        out="metabinner_bins/" + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_metabinner_mags" + "/metabinner_res/metabinner_result.tsv",
    threads:
       long_threads 
    params:
        bin="metabinner_bins/" + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_metabinner_mags",
        curr_dir=CURR_DIR,
        tag=TAG,
    conda:
        envs + "metabinner_env.yaml"
    shell:
        """
        python ../MetaBinner/scripts/gen_kmer.py {input.fasta} 1500 4 
        ../MetaBinner/run_metabinner.sh -a /mnt/disks/MOUNT_DIR/{params.curr_dir}/{input.fasta} -k /mnt/disks/MOUNT_DIR/{params.curr_dir}/{params.tag}{wildcards.sample}/assembly/scaffolds_kmer_4_f1500.csv -o /mnt/disks/MOUNT_DIR/{params.curr_dir}/{params.bin} -d /mnt/disks/MOUNT_DIR/{params.curr_dir}/{input.trim} -t 10 -p /mnt/disks/MOUNT_DIR//MetaBinner/;
        rm -r {params.bin}/metabinner_res/intermediate_result/
        """

rule get_mags:
    input:
        fasta=TAG + "{sample}" + ASS,
        bin="metabinner_bins/" + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_metabinner_mags" + "/metabinner_res/metabinner_result.tsv",
    output:
        folder = directory("metabinner_bins/" + TAG + "{sample}_{mapper}_{sampletype}_metabinner_{seqtype}")
    shell:
        """
        mkdir {output.folder};
        python /mnt/disks/MOUNT_DIR/pipeline/bin_metabinner.py {input.bin} {input.fasta} {output.folder}
        """

rule checkm2_predict_metabinner:
    input:
        folder = "metabinner_bins/" + TAG + "{sample}_{mapper}_{sampletype}_metabinner_{seqtype}"
    output:
        cm2 = directory(checkm2 + TAG + "{sample}_{mapper}_{sampletype}_metabinner_{seqtype}_cm2")
    conda:
        envs + "checkm2_env.yaml"
    threads: short_threads # Adjust the number of threads as necessary
    shell:
        "~/checkm2/bin/checkm2 predict --threads {threads} --input {input.folder} --output-dir {output.cm2} -x fa --force;"
        "rm {output.cm2}/protein_files/*;"
        "rm {output.cm2}/diamond*/*"

    
rule trim_mapping_max:
    input: 
        sorted= covfiles + TAG + "{sample}_"+"{mapper}_all_short.tsv",
    output:
        trim = cc_covfiles +  TAG + "{sample}_"+"{mapper}_all_short_maxbin2.tsv",
    run: 
        import csv

        # Define the input and output files
        input_file = input['sorted']
        output_file = output['trim']

        def filter_rows(input_file, output_file):
            with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
                # Create reader and writer objects
                reader = csv.reader(infile, delimiter='\t')
                writer = csv.writer(outfile, delimiter='\t')

                # Read and write the header
                header = next(reader)
                print(header)
                abridged_header = [header[0]] + [header[i] for i in range(3,len(header),2)]
                writer.writerow(abridged_header)

                # Iterate over the rows
                for row in reader:
                    abridged_row = [row[0]] + [row[i] for i in range(3,len(row),2)]
                    # Split the 'Contig_name' field and extract the length
                    writer.writerow(abridged_row)

        # Call the function with the file paths
        print(input_file, output_file)
        filter_rows(input_file, output_file)


rule maxbin:
    input:
        fasta=TAG + "{sample}" + ASS,
        trim = cc_covfiles +  TAG + "{sample}_"+"{mapper}_all_short_maxbin2.tsv"
    output:
        out=bins + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_maxbin_mags/" + "maxbin_{sample}.summary",
        d = directory(bins + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_maxbin_mags/") 
        #cm2 = directory(checkm2 + TAG + "{sample}_{mapper}_{sampletype}_maxbin_{seqtype}_cm2")
    threads:
        20
    conda:
        envs + "maxbin_env.yaml"
    params:
        bin=bins + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_maxbin_mags",
        curr_dir=CURR_DIR,
        tag=TAG,
    shell:
        """
        run_MaxBin.pl -contig {input.fasta} -out maxbin_{wildcards.sample} -abund {input.trim} -thread 60 -min_contig_length 1500;
        mkdir -p {params.bin};
        mv maxbin_{wildcards.sample}* {params.bin}
        """

rule checkm2_predict_maxbin:
    input:
        folder=bins + TAG + "{sample}_{mapper}_{mapper}_{sampletype}_{seqtype}_maxbin_mags",
    output:
        cm2 = directory(checkm2 + TAG + "{sample}_{mapper}_{sampletype}_maxbin_{seqtype}_cm2")
    conda:
        envs + "checkm2_env.yaml"
    threads: short_threads # Adjust the number of threads as necessary
    shell:
        "~/checkm2/bin/checkm2 predict --threads {threads} --input {input.folder} --output-dir {output.cm2} -x fasta --force;"
        "rm {output.cm2}/protein_files/*;"
        "rm {output.cm2}/diamond*/*"
    
