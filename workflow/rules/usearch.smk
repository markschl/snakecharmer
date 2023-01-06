from os.path import splitext
from lib import file_logging


localrules:
    usearch_stats_paired,


##########################
#### Paired-end
##########################


rule usearch_merge_paired:
    params:
        par=lambda w: cfg[w.name]["settings"]["usearch"]["merge"],
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
            {params.par[program]} {params.par[overlap_ident]} \
            -threads 1 \
            -fastq_maxdiffs {params.par[max_diffs]} \
            2> {log}
        """


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
    log:
        "logs/{name}/usearch/paired/2_trim/{sample}.log",
    conda:
        "envs/cutadapt.yaml"
    group:
        "sample"
    shell:
        """
        $PIPELINE_DIR/workflow/scripts/trim_primers/trim_primers_paired.sh {input.seq} {input.fprimers} {input.rprimers_rev} \
            $(dirname {output[0]}) \
            --error-rate {params.par[max_error_rate]} \
            --overlap {params.par[min_overlap]} \
            --minimum-length {params.minlen} \
             2> {log}
        # create empty zstd archives if file does not exist
        # rename output to contain marker name
        # TODO: kind of complicated procedure
        prefix=processing/{wildcards.name}/usearch/paired/2_trim/{wildcards.sample}/merged
        for marker_comb in {params.primer_comb}; do
            comb=${{marker_comb##*__}}
            if [ ! -f "$prefix/$comb.fastq.zst" ]; then
                echo "No sequences with both forward and reverse primer ($marker_comb) were found in sample {wildcards.sample}" >&2
                echo -n | zstd -cq > "$prefix/$marker_comb.fastq.zst" 2>> {log}
            else
                mv "$prefix/$comb.fastq.zst" "$prefix/$marker_comb.fastq.zst" 2>> {log}
            fi
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
    conda:
        "envs/usearch-vsearch.yaml"
    resources:
        mem_mb=5000,
    shell:
        """
        # This command combines the filtered and de-replicated
        # sample files and will serve as input for the clustering.
        # We therefore de-replicate it again, which is important since
        # UNOISE relies on the size annotations for the -minsize filter
        zstd -dcq {input.good} | 
          vsearch -derep_fulllength - -sizein -sizeout -output - 2> {log} |
          zstd -cq > {output.good}

        # Combine unfiltered (actually, length filtered after trimming)
        # and de-replicated sequences into one file. These will be used
        # for mapping against the denoised sequences in order to create
        # an OTU table. It is important *not* to de-replicate them again here,
        # otherwise many sample labels will be lost, leading to wrong results.
        zstd -dcq {input.all} | zstd -cq > {output.all}
        """


rule usearch_unoise3:
    params:
        par=lambda w: cfg[w.name]["settings"]["usearch"]["unoise"],
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
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    shell:
        """
        if [[ "{params.par[program]}" == "usearch" ]]; then
            zstd -dqf {input}
            f={input}
            usearch -unoise3 ${{f%.zst}} \
            -zotus {output} \
            -minsize {params.par[min_size]} \
            &> {log}
            rm "$f"
        elif [[ "{params.par[program]}" == "vsearch" ]]; then
            zstd -dcqf {input} |
                vsearch -cluster_unoise - \
                    -minsize {params.par[min_size]} \
                    -sizeout \
                    -threads {threads} \
                    -centroids - 2> {log} |
                  vsearch -uchime3_denovo - \
                    -nonchimeras {output} \
                    -sizein \
                    -relabel Zotu \
                    2> {log}
        else
            echo "unknown program: {params.par[program]}"
            exit 1
        fi
        """


rule usearch_make_otutab:
    params:
        par=lambda w: cfg[w.name]["settings"]["usearch"]["otutab"],
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
    resources:
        mem_mb=10000,
        runtime=24 * 60,
    shell:
        """
        # TODO: many arguments, eventually convert into native snakemake Bash/Python script
        $PIPELINE_DIR/workflow/scripts/usearch/make_otutab.sh {params.par[program]} {threads} \
          "{input.uniques}" "{input.denoised}" \
          "{output.tab}" "{output.map}" "{output.bam}" "{output.notmatched}" \
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
        "obtain_taxdb"
    shell:
        """
        $PIPELINE_DIR/workflow/scripts/usearch/convert_utax.sh {input} {output}
        """


rule assign_taxonomy_sintax:
    params:
        par=lambda w: cfg[w.name]["taxonomy"][w.marker][(w.db_name, w.tax_method)],
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
    conda:
        "envs/usearch-vsearch.yaml"
    threads: workflow.cores
    resources:
        mem_mb=5000,
    shell:
        """
        $PIPELINE_DIR/workflow/scripts/usearch/sintax_usearch.sh \
          {params.par[program]} {input.fa} {input.db} {output.tax} \
            -strand both \
            -sintax_cutoff {params.par[confidence]} \
            -threads {threads} \
          &> {log}
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
    run:
        # TODO: turn into script
        from os.path import dirname, basename
        import csv
        from collections import defaultdict


        def read_stats(files, key):
            out = {}
            for f in files:
                with open(f, "r") as h:
                    out[key(f)] = [int(n) for n in next(csv.reader(h, delimiter="\t"))]
            return out


        percent = lambda x, y: round(100 * x / y, 2) if y > 0 else 0.0

        with file_logging(log):
            # read data
            merge = read_stats(input.merge, key=lambda f: basename(dirname(f)))
            trim = read_stats(input.trim, key=lambda f: basename(dirname(dirname(f))))
            flt = read_stats(
                input.filter,
                key=lambda f: (basename(dirname(dirname(f))), basename(dirname(f))),
            )
            flt2 = defaultdict(dict)
            for k, n in flt.items():
                flt2[k[1]][k[0]] = n

            # combine and write to output
            with open(output[0], 'w') as out:
                writer = csv.writer(out, delimiter='\t')
                header = ['sample', 'raw', 'merged', '(%)', 'fwd-primer', '(%)', 'rev-primer', '(%)']
                for p in cfg.primer_combinations_flat:
                    header += [p, "filtered", "(%)"]
                writer.writerow(header)
                for sample in merge:
                    raw, m = merge[sample]
                    m2, tf, tr = trim[sample]
                    assert (m == m2), "Number of sequences in logfiles from read merging and primer trimming does not match: {} vs. {}".format(m, m2)
                    row = [sample, raw, m, percent(m, raw), tf, percent(tf, m), tr, percent(tr, m)]
                    f = flt2[sample]
                    for p in cfg.primer_combinations_flat:
                        kept, removed = f[p]
                        total = kept + removed
                        row += [total, kept, percent(kept, total)]
                    writer.writerow(row)
