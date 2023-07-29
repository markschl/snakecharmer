from os.path import splitext


localrules:
    make_cutadapt_fasta,
    usearch_stats_paired,


##########################
#### Paired-end
##########################


rule usearch_merge_paired:
    params:
        par=lambda w: cfg[w.name]["settings"]["usearch"]["merge"],
        usearch_bin=config['software']['usearch']['binary']
    input:
        expand(
            "input/grouped/paired/{{sample}}/{{sample}}_R{read}.fastq.gz",
            read=[1, 2],
        ),
    output:
        merged="processing/{name}/usearch/paired/1_merged/{sample}/{sample}.fastq.zst",
        r1="processing/{name}/usearch/paired/1_merged/{sample}/{sample}_notmerged_R1.fastq.zst",
        r2="processing/{name}/usearch/paired/1_merged/{sample}/{sample}_notmerged_R2.fastq.zst",
        stats="processing/{name}/usearch/paired/1_merged/{sample}/{sample}_stats.txt",
    log:
        "logs/{name}/usearch/paired/1_merge/{sample}.log",
    conda:
        "envs/usearch-vsearch.yaml"
    group:
        "sample"
    shell:
        """
        out={output.merged}
        $PIPELINE_DIR/workflow/scripts/usearch/merge_paired.sh \
            {input[0]} {input[1]} ${{out%.fastq.zst}} \
            "{params.par[program]}" "{params.usearch_bin}" {params.par[overlap_ident]} \
            -threads 1 \
            -fastq_maxdiffs {params.par[max_diffs]} \
            2> {log}
        """



rule make_cutadapt_fasta:
    input:
        yaml='processing/primers/primers.yaml',
    output:
        forward="processing/primers/forward.fasta",
        reverse_rev="processing/primers/reverse_rev.fasta"
    log:
        "logs/make_primer_fasta.log",
    conda:
        "envs/consensus.yaml"
    script:
        "../scripts/make_primer_fasta.py"


rule trim_primers_paired:
    params:
        par=lambda w: cfg[w.name]["settings"]["primers"]["trim_settings"],
        minlen=lambda w: cfg[w.name]["settings"]["filter"]["min_length"],
        primer_comb=cfg.primer_combinations_flat,
    input:
        fprimers="processing/primers/forward.fasta",
        rprimers_rev="processing/primers/reverse_rev.fasta",
        seq="processing/{name}/usearch/paired/1_merged/{sample}/{sample}.fastq.zst",
    output:
        fwd_log="processing/{name}/usearch/paired/2_trim/{sample}/merged/{sample}_fwd.log",
        rev_log="processing/{name}/usearch/paired/2_trim/{sample}/merged/{sample}_rev.log",
        stats="processing/{name}/usearch/paired/2_trim/{sample}/merged/{sample}_stats.txt",
        by_primers=expand(
            "processing/{{name}}/usearch/paired/2_trim/{{sample}}/merged/{primers}.fastq.zst",
            primers=cfg.primer_combinations_flat,
        ),
        short="processing/{name}/usearch/paired/2_trim/{sample}/merged/too_short.fastq.zst",
    log:
        "logs/{name}/usearch/paired/2_trim/{sample}.log",
    conda:
        "envs/cutadapt.yaml"
    group:
        "sample"
    shell:
        """
        $PIPELINE_DIR/workflow/scripts/trim_primers/trim_primers_paired.sh {input.seq} {input.fprimers} {input.rprimers_rev} \
            $(dirname {output[stats]}) \
            {params.minlen} \
            --error-rate {params.par[max_error_rate]} \
            --overlap {params.par[min_overlap]} \
             2> {log}
        # compress and rename to contain marker name
        # TODO: kind of complicated procedure
        prefix=processing/{wildcards.name}/usearch/paired/2_trim/{wildcards.sample}/merged
        shopt -s nullglob  # TODO: assuming Bash, eventually we may convert to Python script
        for marker_comb in {params.primer_comb}; do
            comb=${{marker_comb##*__}}
            if [ -f "$prefix/$comb.fastq" ]; then
                zstd --rm -qf "$prefix/$comb.fastq" -o "$prefix/$marker_comb.fastq.zst" 2>> {log}
            else
                echo "No sequences with both forward and reverse primer ($marker_comb) were found in sample {wildcards.sample}" >&2
                echo -n | zstd -cq > "$prefix/$marker_comb.fastq.zst" 2>> {log}
            fi
            for f in "$prefix/"no_adapter...*.fastq; do zstd --rm -qf "$f"; done
            for f in "$prefix/"*...no_adapter.fastq; do zstd --rm -qf "$f"; done
        done
        """


rule usearch_filter_derep_paired:
    params:
        maxee=lambda w: cfg[w.name]["settings"]["usearch"]["filter"]["max_error_rate"],
        minlen=lambda w: cfg[w.name]["settings"]["filter"]["min_length"],  # not absolutely required (done after trimming)
    input:
        "processing/{name}/usearch/paired/2_trim/{sample}/merged/{primers}.fastq.zst",
    output:
        all="processing/{name}/usearch/paired/3_filter_derep/{primers}/{sample}/{sample}_all_uniques.fasta.zst",
        good="processing/{name}/usearch/paired/3_filter_derep/{primers}/{sample}/{sample}_good_uniques.fasta.zst",
        stats="processing/{name}/usearch/paired/3_filter_derep/{primers}/{sample}/{sample}_stats.txt",
    log:
        "logs/{name}/usearch/paired/3_filter_derep/{primers}/{sample}.log",
    conda:
        "envs/usearch-vsearch.yaml"
    group:
        "sample"
    resources:
        mem_mb=2000,
    shell:
        """
        # TODO: review posible defaults for fastq_qmax
        all={output.all}
        $PIPELINE_DIR/workflow/scripts/usearch/filter_derep.sh {input} {wildcards.sample} \
            ${{all%_all_uniques.fasta.zst}} \
            -fastq_maxee_rate {params.maxee} \
            -fastq_minlen {params.minlen} \
            -fastq_qmax 42 \
             2> {log}
        """


##########################
#### Denoising/clustering
##########################


rule usearch_collect_derep:
    input:
        good=lambda w: expand(
            "processing/{{name}}/usearch/{{strategy}}/3_filter_derep/{{primers}}/{sample}/{sample}_good_uniques.fasta.zst",
            sample=cfg.sample_names[w.strategy],
        ),
        all=lambda w: expand(
            "processing/{{name}}/usearch/{{strategy}}/3_filter_derep/{{primers}}/{sample}/{sample}_all_uniques.fasta.zst",
            sample=cfg.sample_names[w.strategy],
        ),
    output:
        good="processing/{name}/usearch/{strategy}/4_unique/{primers}/good_uniques.fasta.zst",
        all="processing/{name}/usearch/{strategy}/4_unique/{primers}/all_uniques.fasta.zst",
    log:
        "logs/{name}/usearch/{strategy}/4_collect_derep/{primers}.log",
    group:
        "denoise"
    conda:
        "envs/usearch-vsearch.yaml"
    resources:
        mem_mb=5000,
    shell:
        """
        # This command combines the filtered and de-replicated
        # sample files and will serve as input for the clustering.
        # We therefore de-replicate it again, which is important since
        # UNOISE assumes unique sequneces with size annotations.
        # TODO: the output is sorted by size (see vsearch docs).
        #   -> Make sure that this stays the same when updating to
        #      future versions.
        zstd -dcq {input.good} | 
          vsearch -derep_fulllength - -sizein -sizeout -output - 2> {log} |
          zstd -cq > {output.good}

        # Combine unfiltered (actually, length filtered after trimming)
        # and de-replicated sequences into one file. These will be used
        # for mapping against the denoised sequences to create the OTU table.
        # It is important *not* to de-replicate them again here,
        # otherwise many sample labels will be lost, leading to wrong results.
        zstd -dcq {input.all} | zstd -cq > {output.all}
        """


rule usearch_unoise3:
    params:
        par=lambda w: cfg[w.name]["settings"]["usearch"]["unoise"],
        usearch_bin=config['software']['usearch']['binary']
    input:
        "processing/{name}/usearch/{strategy}/4_unique/{primers}/good_uniques.fasta.zst",
    output:
        "results/{name}/pipeline_usearch_unoise3/{primers}/{strategy}/denoised.fasta",
    log:
        "logs/{name}/usearch/{strategy}/pipeline_usearch_unoise3/{primers}/unoise3.log",
    conda:
        "envs/usearch-vsearch.yaml"
    group:
        "denoise"
    threads:
        # VSEARCH works in parallel (although cores seem to be used only ~50%) while
        # USEARCH v11 does not appear to use more than 1 thread
        # TODO: further validate VSEARCH threads setting
        lambda w: int(workflow.cores*1.5) if cfg[w.name]["settings"]["usearch"]["unoise"]["program"] == "vsearch" else 1
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    shell:
        """
        if [[ "{params.par[program]}" == "usearch" ]]; then
            zstd -dqf {input}
            f={input}
            "{params.usearch_bin}" -unoise3 ${{f%.zst}} \
                -zotus {output} \
                -minsize {params.par[min_size]} \
                &> {log}
            rm "${{f%.zst}}"
        elif [[ "{params.par[program]}" == "vsearch" ]]; then
            # following code from https://github.com/torognes/vsearch/pull/283
            {{
              zstd -dcqf {input} |
                stdbuf -eL vsearch --cluster_unoise - \
                    --minsize {params.par[min_size]} \
                    --sizein \
                    --sizeout \
                    --threads {threads} \
                    --centroids - |
                  stdbuf -eL vsearch --sortbysize - --output - |
                  stdbuf -eL vsearch --uchime3_denovo - \
                    --nonchimeras - \
                    --sizein \
                    --relabel Zotu |
                  st upper --wrap 80 `# convert masked letters to uppercase (alternative: use --qmask none)` \
                  > {output}
            }} 2> {log}
        else
            echo "unknown program: {params.par[program]}"
            exit 1
        fi
        """


rule usearch_make_otutab:
    params:
        par=lambda w: cfg[w.name]["settings"]["usearch"]["otutab"],
        usearch_bin=config['software']['usearch']['binary']
    input:
        denoised="results/{name}/pipeline_usearch_{cluster}/{primers}/{strategy}/denoised.fasta",
        uniques="processing/{name}/usearch/{strategy}/4_unique/{primers}/all_uniques.fasta.zst",
    output:
        tab="results/{name}/pipeline_usearch_{cluster}/{primers}/{strategy}/denoised_otutab.txt.gz",
        map="results/{name}/pipeline_usearch_{cluster}/{primers}/{strategy}/denoised_search.txt.gz",
        notmatched="processing/{name}/usearch/{strategy}/pipeline_usearch_{cluster}/{primers}/make_otutab/notmatched.fasta.zst",
        bam="processing/{name}/usearch/{strategy}/pipeline_usearch_{cluster}/{primers}/make_otutab/mapping.bam",
    threads: workflow.cores
    log:
        "logs/{name}/usearch/{strategy}/pipeline_usearch_{cluster}/{primers}/make_otutab.log",
    conda:
        "envs/vsearch-samtools.yaml"
    group:
        "otutab"
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    shell:
        """
        # TODO: many arguments, eventually convert into native snakemake Bash/Python script
        $PIPELINE_DIR/workflow/scripts/usearch/make_otutab.sh \
          {params.par[program]} "{params.usearch_bin}" \
          {threads} \
          "{input.uniques}" "{input.denoised}" \
          "{output.tab}" "{output.map}" \
          "{output.bam}" "{output.notmatched}" \
         -id {params.par[ident_threshold]} \
         -maxaccepts {params.par[maxaccepts]} \
         -maxrejects {params.par[maxrejects]} \
         2> {log}
        """


##########################
#### Taxonomy
##########################


rule convert_taxdb_utax:
    input:
        "refdb/taxonomy/{db}/filtered/{defined}.fasta.zst",
    output:
        "refdb/taxonomy/{db}/formatted/utax/{defined}.fasta.zst",
    conda:
        "envs/basic.yaml"
    group:
        "taxonomy"
    shell:
        """
        $PIPELINE_DIR/workflow/scripts/usearch/convert_utax.sh {input} {output}
        """


rule assign_taxonomy_sintax:
    params:
        par=lambda w: cfg[w.name]["taxonomy"][w.marker][(w.db_name, w.tax_method)],
        usearch_bin=config['software']['usearch']['binary']
    input:
        fa="results/{name}/{pipeline}/{marker}__{primers}/{strategy}/denoised.fasta",
        db=lambda w: "refdb/taxonomy/{{marker}}/{db_name}/formatted/utax/{defined}.fasta.zst".format(
            **cfg[w.name]["taxonomy"][w.marker][(w.db_name, w.tax_method)]
        ),
    output:
        tax="results/{name}/{pipeline}/{marker}__{primers}/{strategy}/taxonomy/{db_name}-sintax_usearch-{tax_method}.txt.gz",
        sintax="results/{name}/{pipeline}/{marker}__{primers}/{strategy}/taxonomy/sintax/{db_name}-sintax_usearch-{tax_method}.txt.gz",
    log:
        "logs/{name}/other/{strategy}/{pipeline}/{marker}__{primers}/taxonomy_sintax/{db_name}-{tax_method}.log",
    group:
        "taxonomy"
    conda:
        "envs/usearch-vsearch.yaml"
    threads:
        # VSEARCH works in parallel (although cores seem to be used only ~50%) while
        # USEARCH v11 does not appear to use more than 1 thread
        lambda w: workflow.cores if cfg[w.name]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["program"] == "vsearch" else 1,
    resources:
        mem_mb=5000,
    shell:
        """
        bin="{params.par[program]}"
        if [ $bin = "usearch" ]; then
            bin="{params.usearch_bin}"
        fi
        $PIPELINE_DIR/workflow/scripts/usearch/sintax_usearch.sh \
          {params.par[program]} "$bin" {input.fa} {input.db} {output.tax} \
            -strand both \
            -sintax_cutoff {params.par[confidence]} \
            -threads {threads} \
          2> {log}
        """


##########################
#### QC
##########################


rule usearch_multiqc_paired:
    input:
        fastqc=rules.multiqc_fastqc.input,
        cutadapt=expand(
            "processing/{{name}}/usearch/paired/2_trim/{sample}/merged/{sample}_{dir}.log",
            sample=cfg.sample_names["paired"],
            dir=["fwd", "rev"],
        ),
    output:
        "results/{name}/pipeline_usearch_{cluster}/_validation/multiqc/multiqc_report.html",
    log:
        "logs/{name}/usearch/paired/pipeline_usearch_{cluster}/multiqc_paired.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        cutadapt_dir=$(dirname $(dirname $(dirname {input.cutadapt[0]})))
        multiqc -f -m fastqc -m cutadapt -o $(dirname {output}) results/_validation/fastqc $cutadapt_dir 2> {log}
        """


# TODO: parse both paired and single into one sample_report.tsv
rule usearch_stats_paired:
    params:
        primer_combinations=cfg.primer_combinations_flat
    input:
        merge=expand(
            "processing/{{name}}/usearch/paired/1_merged/{sample}/{sample}_stats.txt",
            sample=cfg.sample_names["paired"],
        ),
        trim=expand(
            "processing/{{name}}/usearch/paired/2_trim/{sample}/merged/{sample}_stats.txt",
            sample=cfg.sample_names["paired"],
        ),
        filter=expand(
            "processing/{{name}}/usearch/paired/3_filter_derep/{primers}/{sample}/{sample}_stats.txt",
            primers=cfg.primer_combinations_flat,
            sample=cfg.sample_names["paired"],
        ),
    output:
        "results/{name}/pipeline_usearch_{cluster}/_validation/sample_report.tsv",
    log:
        "logs/{name}/usearch/paired/pipeline_usearch_{cluster}/sample_report.log",
    script:
        "../scripts/usearch/stats_paired.py"
