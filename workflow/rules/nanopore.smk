samples_nano = get_samples('nanopore')

use rule rename_ref as rename_ref_nano with:
    input: 
        config["ref"]["nanopore"]
    output: "nanopore/ref.fasta"
    log: "logs/nanopore/rename_ref.log"

# create abundance matrix with minimap
rule minimap2_index:
    input: rules.rename_ref_nano.output,
    output: temp("nanopore/ref.mmi")
    params:
        index_size = "4G",
    conda: "../envs/minimap2.yaml"
    log: "logs/nanopore/mini2index.log"
    shell: "minimap2 -I {params.index_size} -d {output} {input} 2> {log}"

rule samtools_dict:
    input: rules.rename_ref_nano.output,
    output: temp("nanopore/ref.dict")
    conda: "../envs/samtools.yaml"
    log: "logs/nanopore/samdict.log"
    shell: "samtools dict {input} | cut -f1-3 > {output} 2> {log}"

rule minimap2:
    input:
        fq = lambda wc: get_fqs(wc)[0],
        mmi = rules.minimap2_index.output,
        dict = rules.samtools_dict.output,
    output: temp("nanopore/mapped/{sample}.sam")
    params:
        x = "map-ont",
    conda: "../envs/minimap2.yaml"
    log: "logs/nanopore/minimap2/{sample}.log"
    threads: config["threads"]["large"]
    shell:
        """
        minimap2 -t {threads} -ax {params.x} --secondary=no {input.mmi} {input.fq} 2> {log} | \
        grep -v "^@" | cat {input.dict} - > {output} 2>> {log}
        """

rule samtools_bamsort:
    input: rules.minimap2.output
    output: temp("nanopore/mapped/{sample}.sorted.bam")
    params:
        prefix = "nanopore/mapped/tmp.{sample}",
        m = "3G",
    conda: "../envs/samtools.yaml"
    log: "logs/nanopore/samsort/{sample}.log"
    shell:
        "samtools view -F 3584 -b {input} | samtools sort - -T {params.prefix} -m {params.m} -o {output} 2>{log}"

use rule samtools_index as samtools_index_nano with:        
    input: rules.samtools_bamsort.output
    output: temp("nanopore/mapped/{sample}.sorted.bam.bai")
    log: "logs/nanopore/samindex/{sample}.log"

# biom format header
use rule ids as ids_nano with:
    input:
        bam = expand("nanopore/mapped/{sample}.sorted.bam", sample=samples_nano),
        bai = expand("nanopore/mapped/{sample}.sorted.bam.bai", sample=samples_nano)
    output: 
        temp("nanopore/ids")
    log: 
        "logs/nanopore/ids_nano.log"
       
use rule counts as counts_nano with:
    input:
        bam = "nanopore/mapped/{sample}.sorted.bam",
        bai = "nanopore/mapped/{sample}.sorted.bam.bai"
    output: 
        temp("nanopore/mapped/{sample}.count")
    log: 
        "logs/nanopore/counts/{sample}.log"

use rule matrix as matrix_nano with:
    input:
        ids = rules.ids_nano.output,
        counts = expand("nanopore/mapped/{sample}.count", sample=samples_nano),
    output: 
        "matrix_nanopore.tsv"
