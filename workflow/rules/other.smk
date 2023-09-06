

localrules:
    clean_itsx,
    clean_cmp,


rule itsx:
    params:
        par=config["ITSx"],
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}_{layout}/ITS__{primers}/denoised.fasta",
    output:
        pos="results/{workflow}/workflow_{cluster}/{run}_{layout}/ITS__{primers}/ITSx/out.positions.txt",
    log:
        "logs/results/{workflow}/workflow_{cluster}/{run}_{layout}/ITS__{primers}/ITSx.log",
    group:
        "ITS"
    conda:
        "envs/itsx.yaml"
    threads: min(workflow.cores, 8)  # too many does not make sense
    shell:
        """
        pos="{output.pos}"
        ITSx -i {input} -o "${{pos%.positions.txt}}" \
          -t {params.par[organism_groups]} \
          -E {params.par[e-value]} \
          --allow_single_domain {params.par[allow_single_domain]} \
          --complement {params.par[complement]} \
          --heuristics {params.par[heuristics]} \
          --graphical {params.par[graphical]} \
          --fasta {params.par[fasta]} \
          --preserve {params.par[preserve]} \
          --save_regions {params.par[save_regions]} \
          --partial {params.par[partial]} \
          --cpu {threads} \
          2> {log}
        """


rule vsearch_global:
    params:
        db=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "file"),
        threshold=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "ident_threshold"),
        maxaccepts=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "maxaccepts"),
        maxrejects=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "maxrejects"),
        maxhits=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "maxhits"),
    input:
        otus="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/denoised.fasta",
    output:
        map="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}.txt",
        bam="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}.bam",
        denoised_notmatched="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}_denoised_notmatched.fasta.gz",
        notmatched="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}_notmatched.fasta.gz",
    log:
        "logs/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp_{comparison}.log",
    group:
        "cmp"
    # TODO: resources? usually pretty fast
    threads: max(1, int(workflow.cores / 2))
    conda:
        "envs/vsearch-samtools.yaml"
    shell:
        """
        bam={output.bam}
        sam=${{bam%.*}}
        notmatched={output.notmatched}; notmatched=${{notmatched%.gz}}
        denoised_notmatched={output.denoised_notmatched}; denoised_notmatched=${{denoised_notmatched%.gz}}
        vsearch -usearch_global {input.otus} -db "{params.db}" \
            -userout {output.map} \
            -samout $sam \
            -userfields 'query+target+id' \
            -notmatched $denoised_notmatched \
            -dbnotmatched $notmatched \
            -threads {threads} \
            -maxaccepts {params.maxaccepts} \
            -maxrejects {params.maxrejects} \
            -maxhits {params.maxhits} \
            -id {params.threshold} &> {log}
        # compress not-matched files
        gzip -nf $notmatched $denoised_notmatched 2> {log}
        # make BAM file
        rm -f "{params.db}.fai" $bam.bai
        samtools view -T "{params.db}" -b $sam |
          samtools sort -@ {threads} > $bam 2> {log}
        rm -f $sam "{params.db}.fai"
        samtools index $bam 2> {log}
        """


rule clean_itsx:
    shell:
        "rm -Rf results/*/workflow_*/*/*/ITSx"


rule clean_cmp:
    shell:
        "rm -Rf results/*/workflow_*/*/*/cmp"
