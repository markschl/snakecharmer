
localrules:
    clean_cmp,


rule vsearch_global:
    params:
        db=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "file"),
        threshold=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "ident_threshold"),
        maxaccepts=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "maxaccepts"),
        maxrejects=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "maxrejects"),
        maxhits=lambda w: with_default(cfg[w.workflow]["settings"]["compare"], w.comparison, "maxhits"),
    input:
        otus="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/clusters_fwd.fasta",
    output:
        map="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}.txt",
        bam="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}.bam",
        clusters_notmatched="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}_clusters_notmatched.fasta.gz",
        notmatched="results/{workflow}/workflow_{cluster}/{run}_{layout}/{primers}/cmp/{comparison}_notmatched.fasta.gz",
    log:
        "logs/{workflow}/{run}_{layout}/{primers}/{cluster}_cmp_{comparison}.log",
    group:
        "cmp"
    # TODO: resources? usually pretty fast
    threads: max(1, int(workflow.cores / 2))
    conda:
        "envs/vsearch-samtools.yaml"
    shell:
        """
        exec &> {log}
        set -xeuo pipefail
        bam="{output.bam}"
        sam="${{bam%.*}}.sam"
        notmatched="{output.notmatched}"
        notmatched="${{notmatched%.gz}}"
        clusters_notmatched="{output.clusters_notmatched}"
        clusters_notmatched="${{clusters_notmatched%.gz}}"
        vsearch -usearch_global "{input.otus}" -db "{params.db}" \
            -userout "{output.map}" \
            -samout "$sam" \
            -userfields 'query+target+id' \
            -notmatched "$clusters_notmatched" \
            -dbnotmatched "$notmatched" \
            -threads {threads} \
            -maxaccepts {params.maxaccepts} \
            -maxrejects {params.maxrejects} \
            -maxhits {params.maxhits} \
            -id {params.threshold}
        # compress not-matched files
        gzip -nf "$notmatched" "$clusters_notmatched"
        if [ -s "$sam" ]; then
            # make BAM file
            rm -f "{params.db}.fai" "$bam.bai"
            samtools view -T "{params.db}" -b "$sam" |
            samtools sort -@ {threads} > "$bam"
            samtools index "$bam"
        else
            echo -n > "$bam"
        fi
        rm -f "$sam" "{params.db}.fai"
        """


rule clean_cmp:
    shell:
        "rm -Rf results/*/workflow_*/*/*/cmp"
