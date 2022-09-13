samples_illu = get_samples('illumina')

rule rename_ref:
    input: config["ref"]["illumina"]
    output: "illumina/ref.fasta"
    log: "logs/illumina/rename_ref.log"
    run:
        with open(output[0], "w") as out:
            with open (input[0], "r") as inp:
                i = 1
                for line in inp:
                    if line.startswith(">"):
                        line = ">id_" + str(i) + "\n"
                        i += 1 
                    out.write(line)

rule bwt2_build:
    input: 
        ref = rules.rename_ref.output,
    output:
        multiext("illumina/ref", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: "logs/illumina/bwt2build.log"
    params:
        extra=""  # optional parameters
    threads: config["threads"]["normal"]
    wrapper: "v1.12.2/bio/bowtie2/build"

rule bwt2_map:
    input:
        idx=rules.bwt2_build.output,
        sample = lambda wc: get_fqs(wc),
    output: temp("illumina/{sample}.bam")
    log: "logs/illumina/bwt2map/{sample}.log"
    params:
        index=os.getcwd() + "/illumina/ref",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: config["threads"]["normal"]  # Use at least two threads
    wrapper: "v1.12.2/bio/bowtie2/align"

rule samtools_sort:
    input: rules.bwt2_map.output
    output: temp("illumina/{sample}.sorted.bam")
    conda: "../envs/samtools.yaml"
    log: "logs/illumina/samsort/{sample}.log"
    shell:
        "(samtools view -h {input} | {workflow.basedir}/scripts/read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o {output} -) 2>> {log}"

rule samtools_index:        
    input: rules.samtools_sort.output
    output: temp("illumina/{sample}.sorted.bam.bai")
    log: "logs/illumina/samindex/{sample}.log"
    threads:  # Samtools takes additional threads through its option -@
        config["threads"]["normal"]     # This value - 1 will be sent to -@
    wrapper: "v1.12.2/bio/samtools/index"

rule ids:
    input:
        bam=expand("illumina/{sample}.sorted.bam", sample=samples_illu),
        bai=expand("illumina/{sample}.sorted.bam.bai", sample=samples_illu)
    output: temp("illumina/ids")
    conda: "../envs/samtools.yaml"
    log: "logs/illumina/ids.log"
    shell:
        """
        echo '#OTU ID' > {output}
        samtools idxstats $(ls {input.bam} | head -n 1) | grep -v "*" | cut -f1 >> {output}
        """
       
rule counts:
    input:
        bam="illumina/{sample}.sorted.bam",
        bai="illumina/{sample}.sorted.bam.bai"
    output: temp("illumina/counts/{sample}.count")
    conda: "../envs/samtools.yaml"
    log: "logs/illumina/counts/{sample}.log"
    shell:
        """
        echo '{wildcards.sample}' > {output}
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 >> {output}
        """

rule matrix:
    input:
        ids=rules.ids.output,
        counts=expand("illumina/counts/{sample}.count", sample=samples_illu),
    output: "matrix_illumina.tsv"
    log: "logs/illumina/count_matrix.log"
    shell:
        """
        cp {input.ids} {output}
        # split if exceed 1000 (max 1024 by default in most OS)
        while read fs; do
            paste {output} $fs > {output}.tmp
            # redirection first
            mv {output}.tmp {output}
        done < <(echo {input.counts} | xargs -n 1000)
        """
