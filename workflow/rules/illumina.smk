samples_illu = get_samples(seq_type='illumina')

rule bwt2_build:
    input: config["ref"]["illumina"]
    output:
        multiext("illumina/ref", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: "logs/illumina/bwt2build.log"
    conda: "../envs/bwt2sam.yaml" # libtbb.so.2 incompetence due to conda-forge updates
    params:
        extra=""  # optional parameters
    threads: THREADS
    wrapper:
        "v0.75.0/bio/bowtie2/build"

rule bwt2_map:
    input:
        index=rules.bwt2_build.output,
        sample=get_fqs,
    output: temp("illumina/{sample}.bam")
    log:
        "logs/illumina/bwt2map/{sample}.log"
    conda:
        "../envs/bwt2sam.yaml" # libtbb.so.2 incompetence due to conda-forge updates
    params:
        index=os.get_cwd() + "/illumina/ref",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: THREADS  # Use at least two threads
    wrapper:
        "v0.75.0/bio/bowtie2/align"

rule samtools_sort:
    input: rules.bwt2_map.output
    output: temp("illumina/{sample}.sorted.bam")
    conda: "../envs/bwt2sam.yaml"
    log: "logs/illumina/samsort/{sample}.log"
    shell:
        "(samtools view -h {input} | scripts/read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o {output} -) 2>> {log}"

rule samtools_index:        
    input: rules.samtools_sort.output
    output: temp("illumina/{sample}.sorted.bam.bai")
    log: "logs/illumina/samindex/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        THREADS     # This value - 1 will be sent to -@
    wrapper:
        "v0.75.0/bio/samtools/index"

rule ids_illu:
    input:
        bam=expand("illumina/{sample}.sorted.bam", sample=samples_illu),
        bai=expand("illumina/{sample}.sorted.bam.bai", sample=samples_illu)
    output:
        temp("illumina/ids")
    conda:
        "../envs/bwt2sam.yaml"
    log:
        "logs/illumina/ids.log"
    shell:
        """
        echo 'motus' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule counts_illu:
    input:
        bam="illumina/{sample}.sorted.bam",
        bai="illumina/{sample}.sorted.bam.bai"
    output: temp("illumina/counts/{sample}.count")
    conda: "../envs/bwt2sam.yaml"
    log: "logs/illumina/counts/{sample}.log"
    shell:
        """
        echo '{wildcards.sample}' > {output}
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 >> {output}
        """

rule matrix_illu:
    input:
        ids=rules.ids_illu.output,
        counts=expand("illumina/counts/{sample}.count", sample=samples_illu),
    output: "matrix_illumina.tsv"
    log: "logs/illumina/count_matrix.log"
    shell:
        """
        paste {input.ids} {input.counts} > {output}
        """