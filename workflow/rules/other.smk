

rule itsx:
    params:
        par=config["ITSx"],
    input:
        fa="results/{name}/{pipeline}/{primers}/{strategy}/denoised.fasta",
    output:
        pos="results/{name}/{pipeline}/{primers}/{strategy}/ITSx/out.positions.txt",
    log:
        "logs/{name}/other/{strategy}/{pipeline}/{primers}/ITSx.log",
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
        par=lambda w: cfg.cmp_files[w.db],
        maxhits=lambda w: cfg.cmp_files[w.db].get("maxhits", 0),  # 0 = unlimited
    input:
        otus="results/{name}/{pipeline}/{primers}/{strategy}/denoised.fasta",
    output:
        map="results/{name}/{pipeline}/{primers}/{strategy}/cmp/{db}.txt",
        bam="results/{name}/{pipeline}/{primers}/{strategy}/cmp/{db}.bam",
        denoised_notmatched="results/{name}/{pipeline}/{primers}/{strategy}/cmp/{db}_denoised_notmatched.fasta.gz",
        notmatched="results/{name}/{pipeline}/{primers}/{strategy}/cmp/{db}_notmatched.fasta.gz",
    log:
        "logs/{name}/other/{strategy}/{pipeline}/{primers}/search_{db}.log",
    # TODO: resources? usually pretty fast
    threads: round(workflow.cores / 2)
    conda:
        "envs/vsearch-samtools.yaml"
    shell:
        """
        bam={output.bam}
        sam=${{bam%.*}}
        db="{params.par[file]}"
        notmatched={output.notmatched}; notmatched=${{notmatched%.gz}}
        denoised_notmatched={output.denoised_notmatched}; denoised_notmatched=${{denoised_notmatched%.gz}}
        vsearch -usearch_global {input.otus} -db "$db" \
            -userout {output.map} \
            -samout $sam \
            -userfields 'query+target+id' \
            -notmatched $denoised_notmatched \
            -dbnotmatched $notmatched \
            -threads {threads} \
            -maxaccepts {params.par[maxaccepts]} \
            -maxrejects {params.par[maxaccepts]} \
            -maxhits {params.maxhits} \
            -id {params.par[ident_threshold]} &> {log}
        # compress not-matched files
        gzip -f $notmatched $denoised_notmatched 2> {log}
        # make BAM file
        rm -f "$db.fai" $bam.bai
        samtools view -T "$db" -b $sam |
          samtools sort -@ {threads} > $bam 2> {log}
        rm -f $sam "$db.fai"
        samtools index $bam 2> {log}
        """
