
localrules:
    make_cutadapt_fasta,
    usearch_stats_paired,


##########################
#### Paired-end
##########################


rule usearch_merge_paired:
    params:
        par=lambda w: cfg[w.workflow]["settings"]["usearch"]["merge"],
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        r1="processing/{workflow}/input/illumina/paired/{run}/{sample}_R1.fastq.gz",
        r2="processing/{workflow}/input/illumina/paired/{run}/{sample}_R2.fastq.gz",
    output:
        merged="processing/{workflow}/usearch/{run}/paired/1_merged/{sample}/{sample}.fastq.zst",
        r1="processing/{workflow}/usearch/{run}/paired/1_merged/{sample}/{sample}_notmerged_R1.fastq.zst",
        r2="processing/{workflow}/usearch/{run}/paired/1_merged/{sample}/{sample}_notmerged_R2.fastq.zst",
        stats="processing/{workflow}/usearch/{run}/paired/1_merged/{sample}/{sample}_stats.txt",
    log:
        "logs/{workflow}/usearch/{run}/paired/1_merge/{sample}.log",
    conda:
        "envs/usearch-vsearch.yaml"
    group:
        "sample"
    shell:
        """
        out={output.merged}
        $PIPELINE_DIR/workflow/scripts/usearch/merge_paired.sh \
            {input.r1} {input.r2} ${{out%.fastq.zst}} \
            "{params.par[program]}" "{params.usearch_bin}" {params.par[overlap_ident]} \
            -threads 1 \
            -fastq_maxdiffs {params.par[max_diffs]} \
            2> {log}
        """


rule make_cutadapt_fasta:
    input:
        yaml="processing/primers/primers.yaml",
    output:
        forward="processing/primers/forward.fasta",
        reverse_rev="processing/primers/reverse_rev.fasta",
    log:
        "logs/make_primer_fasta.log",
    conda:
        "envs/consensus.yaml"
    script:
        "../scripts/make_primer_fasta.py"


rule trim_primers_paired:
    params:
        par=lambda w: cfg[w.workflow]["settings"]["primers"]["trim_settings"],
        minlen=lambda w: cfg[w.workflow]["settings"]["filter"]["min_length"],
        primer_comb=cfg.primer_combinations_flat,
    input:
        fprimers="processing/primers/forward.fasta",
        rprimers_rev="processing/primers/reverse_rev.fasta",
        seq="processing/{workflow}/usearch/{run}/paired/1_merged/{sample}/{sample}.fastq.zst",
    output:
        fwd_log="processing/{workflow}/usearch/{run}/paired/2_trim/{sample}/merged/{sample}_fwd.log",
        rev_log="processing/{workflow}/usearch/{run}/paired/2_trim/{sample}/merged/{sample}_rev.log",
        stats="processing/{workflow}/usearch/{run}/paired/2_trim/{sample}/merged/{sample}_stats.txt",
        by_primers=expand(
            "processing/{{workflow}}/usearch/{{run}}/paired/2_trim/{{sample}}/merged/{primers}.fastq.zst",
            primers=cfg.primer_combinations_flat,
        ),
        short="processing/{workflow}/usearch/{run}/paired/2_trim/{sample}/merged/too_short.fastq.zst",
    log:
        "logs/{workflow}/usearch/{run}/paired/2_trim/{sample}.log",
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
        # TODO: default to --no-indels? \
        # compress and rename to contain marker name
        # TODO: kind of complicated procedure
        prefix=$(dirname "{output.stats}")
        shopt -s nullglob  # TODO: assuming Bash, eventually we may convert to Python script
        for marker_comb in {params.primer_comb}; do
            comb=${{marker_comb##*__}}
            if [ -f "$prefix/$comb.fastq" ]; then
                zstd --rm -qf "$prefix/$comb.fastq" -o "$prefix/$marker_comb.fastq.zst" 2>> {log}
            else
                echo "No sequences with both forward and reverse primer ($marker_comb) were found in sample {wildcards.sample}" > {log}
                echo -n | zstd -cq > "$prefix/$marker_comb.fastq.zst" 2>> {log}
            fi
        done
        for f in "$prefix/"no_adapter...*.fastq; do zstd --rm -qf "$f"; done
        for f in "$prefix/"*...no_adapter.fastq; do zstd --rm -qf "$f"; done
        """


rule usearch_filter_derep_paired:
    params:
        maxee=lambda w: cfg[w.workflow]["settings"]["usearch"]["filter"]["max_error_rate"],
        minlen=lambda w: cfg[w.workflow]["settings"]["filter"]["min_length"],  # not absolutely required (done after trimming)
    input:
        "processing/{workflow}/usearch/{run}/paired/2_trim/{sample}/merged/{primers}.fastq.zst",
    output:
        all="processing/{workflow}/usearch/{run}/paired/3_filter_derep/{primers}/{sample}/{sample}_all_uniques.fasta.zst",
        good="processing/{workflow}/usearch/{run}/paired/3_filter_derep/{primers}/{sample}/{sample}_good_uniques.fasta.zst",
        stats="processing/{workflow}/usearch/{run}/paired/3_filter_derep/{primers}/{sample}/{sample}_stats.txt",
    log:
        "logs/{workflow}/usearch/{run}/paired/3_filter_derep/{primers}/{sample}.log",
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
        # The 'expand_samples' method depends on the 'get_sample_config' and 'link_samples' checkpoints
        good=lambda wildcards: expand_samples(
            "processing/{workflow}/usearch/{run}/{layout}/3_filter_derep/{primers}/{sample}/{sample}_good_uniques.fasta.zst",
            technology="illumina",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards
        ),
        all=lambda wildcards: expand_samples(
            "processing/{workflow}/usearch/{run}/{layout}/3_filter_derep/{primers}/{sample}/{sample}_all_uniques.fasta.zst",
            technology="illumina",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards
        ),
    output:
        good="processing/{workflow}/usearch/{run}/{layout}/4_unique/{primers}/good_uniques.fasta.zst",
        all="processing/{workflow}/usearch/{run}/{layout}/4_unique/{primers}/all_uniques.fasta.zst",
    log:
        "logs/{workflow}/usearch/{run}/{layout}/4_unique/{primers}.log",
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
        par=lambda w: cfg[w.workflow]["settings"]["usearch"]["unoise"],
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        "processing/{workflow}/usearch/{run}/{layout}/4_unique/{primers}/good_uniques.fasta.zst",
    output:
        "results/{workflow}/workflow_usearch_unoise3/{run}_{layout}/{primers}/denoised.fasta",
    log:
        "logs/{workflow}/usearch/{run}_{layout}/workflow_usearch_unoise3/{primers}/unoise3.log",
    conda:
        "envs/usearch-vsearch.yaml"
    group:
        "denoise"
    # threads:
    # VSEARCH works in parallel (although cores seem to be used only ~50%) while
    # USEARCH v11 does not appear to use more than 1 thread
    # TODO: further validate VSEARCH threads setting
    threads:
        lambda w: int(workflow.cores * 1.5) if cfg[w.workflow]["settings"]["usearch"]["unoise"]["program"] == "vsearch" else 1
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
        par=lambda w: cfg[w.workflow]["settings"]["usearch"]["otutab"],
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        denoised="results/{workflow}/workflow_usearch_{cluster}/{run}_{layout}/{primers}/denoised.fasta",
        uniques="processing/{workflow}/usearch/{run}/{layout}/4_unique/{primers}/all_uniques.fasta.zst",
    output:
        tab="results/{workflow}/workflow_usearch_{cluster}/{run}_{layout}/{primers}/denoised_otutab.txt.gz",
        map="results/{workflow}/workflow_usearch_{cluster}/{run}_{layout}/{primers}/denoised_search.txt.gz",
        notmatched="processing/{workflow}/usearch/{run}/{layout}/workflow_usearch_{cluster}/{primers}/make_otutab/notmatched.fasta.zst",
        bam="processing/{workflow}/usearch/{run}/{layout}/workflow_usearch_{cluster}/{primers}/make_otutab/mapping.bam",
    threads: workflow.cores
    log:
        "logs/{workflow}/usearch/{run}_{layout}/workflow_usearch_{cluster}/{primers}/make_otutab.log",
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
        db="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/qiime.fasta.zst"
    output:
        db="refdb/taxonomy/db_regular_{source_id}/flt_{filter_id}/utax.fasta.zst"
    log:
        "logs/taxdb/filter/db_regular_{source_id}/flt_{filter_id}/convert_utax.log",
    wildcard_constraints:
        source_id = "\w+",
        filter_id = "\w+",
    conda:
        "envs/taxonomy.yaml"
    group:
        "taxonomy"
    script:
        "../scripts/taxonomy/convert_utax.py"


rule assign_taxonomy_sintax:
    params:
        par=lambda w: cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["assign"],
        usearch_bin=config["software"]["usearch"]["binary"],
    input:
        fa="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/denoised.fasta",
        db=lambda w: "refdb/taxonomy/db_{source[preformatted]}_{source[source_id]}/flt_{filter_id}/utax.fasta.zst".format(
            **cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]
        ),
    output:
        tax="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy/{db_name}-sintax_usearch-{tax_method}.txt.gz",
        sintax="results/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy/sintax/{db_name}-sintax_usearch-{tax_method}.txt.gz",
    log:
        "logs/{workflow}/workflow_{cluster}/{run}_{layout}/{marker}__{primers}/taxonomy_sintax/{db_name}-{tax_method}.log",
    group:
        "taxonomy"
    conda:
        "envs/usearch-vsearch.yaml"
    # threads:
    # VSEARCH works in parallel (although cores seem to be used only ~50%) while
    # USEARCH v11 does not appear to use more than 1 thread
    threads:
        lambda w: workflow.cores \
            if cfg[w.workflow]["taxonomy"][w.marker][(w.db_name, w.tax_method)]["assign"]["program"] == "vsearch" \
            else 1
    resources:
        mem_mb=5000,
    shell:
        """
        program="{params.par[program]}"
        usearch_bin="$program"
        if [ "$program" = "usearch" ]; then
            usearch_bin="{params.usearch_bin}"
        fi
        $PIPELINE_DIR/workflow/scripts/usearch/sintax_usearch.sh \
          "$program" "$usearch_bin" {input.fa} {input.db} {output.tax} \
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
        fastqc=lambda wildcards: expand_samples(
            path="input/{technology}/paired/{demux_method}/{run}/fastqc/{sample}_R{read}_fastqc.html",
            demux_method="demux",
            **wildcards
        ),  #demux_method=cfg[w.workflow]["demux_method"],
        cutadapt=lambda wildcards: expand_samples(
            "processing/{workflow}/usearch/{run}/{layout}/2_trim/{sample}/merged/{sample}_{dir}.log",
            technology="illumina",
            dir=["fwd", "rev"],
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards,
        ),
    output:
        "results/{workflow}/validation/multiqc_usearch/multiqc_report.html",
    log:
        "logs/{workflow}/multiqc_usearch.log",
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        indir=input
        outdir="$(dirname {output})"
        cutadapt_dir=$(dirname $(dirname $(dirname {input.cutadapt[0]})))
        echo $indir $outdir $cutadapt_dir
        multiqc -f -m fastqc -m cutadapt -o "$outdir" "$indir" "$cutadapt_dir" 2> {log}
        """


# TODO: parse both paired and single into one sample_report.tsv
rule usearch_stats_paired:
    params:
        primer_combinations=cfg.primer_combinations_flat,
    input:
        merge=lambda wildcards: expand_samples(
            "processing/{{workflow}}/usearch/{run}/{layout}/1_merged/{sample}/{sample}_stats.txt",
            technology="illumina",
            layout="paired",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards,
        ),
        trim=lambda wildcards: expand_samples(
            "processing/{{workflow}}/usearch/{run}/{layout}/2_trim/{sample}/merged/{sample}_stats.txt",
            technology="illumina",
            layout="paired",
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards,
        ),
        filter=lambda wildcards: expand_samples(
            "processing/{{workflow}}/usearch/{run}/{layout}/3_filter_derep/{primers}/{sample}/{sample}_stats.txt",
            layout="paired",
            technology="illumina",
            primers=cfg.primer_combinations_flat,
            pool=cfg[wildcards.workflow]["settings"]["pool_raw"],
            **wildcards,
        ),
    output:
        "results/{workflow}/workflow_usearch_{cluster}/{run}_paired/sample_report.tsv",
    log:
        "logs/{workflow}/workflow_usearch_{cluster}/{run}_paired/sample_report.log",
    script:
        "../scripts/usearch/stats_paired.py"
