
# Possible taxonomic ranks with their corresponding QIIME prefix,
# following https://github.com/bokulich-lab/RESCRIPt
RANK_TRANS <- c(
    'domain' = 'd',
    'superkingdom' = 'sk',
    'kingdom' = 'k',
    'subkingdom' = 'ks',
    'superphylum' = 'sp',
    'phylum' = 'p',
    'subphylum' = 'ps',
    'infraphylum' = 'pi',
    'superclass' = 'sc',
    'class' = 'c',
    'subclass' = 'cs',
    'infraclass' = 'ci',
    'cohort' = 'co',
    'superorder' = 'so',
    'order' = 'o',
    'suborder' = 'os',
    'infraorder' = 'oi',
    'parvorder' = 'op',
    'superfamily' = 'sf',
    'family' = 'f',
    'subfamily' = 'fs',
    'tribe' = 't',
    'subtribe' = 'ts',
    'genus' = 'g',
    'subgenus' = 'gs',
    'species group' = 'ss',
    'species subgroup' = 'sgs',
    'species' = 's',
    'subspecies' = 'ssb',
    'forma' = 'for'
)

make_lineage <- function(ids, total_ranks) {
  artificial_ranks = NULL
  if (is.null(ids[[1]]$rank)) {
    artificial_ranks <- paste0("r", 0:(total_ranks-1))
  }
  sapply(seq_along(ids), function(i) {
        x = ids[[i]]
        n_ranks = length(x$taxon)
        lineage = x$taxon[2:n_ranks]
        if (!is.null(artificial_ranks)) {
          ranks_short = artificial_ranks[2:n_ranks]
        } else {
          ranks = x$rank[2:n_ranks]
          ranks_short <- unname(RANK_TRANS[ranks])
          is_na <- is.na(ranks_short)
          ranks_short[is_na] <- ranks[is_na]
        }
        setNames(
            x$confidence[n_ranks],
            paste(paste(ranks_short, lineage, sep="__"), collapse=";")
        )
    })
}

get_default <- function(list, item, default=NULL) {
    value <- list[[item]]
    if (is.null(value)) default else value
}

idtaxa_assign <- function(seqfile, db, tax_out, threshold, processors=1, rand_seed=NULL, ...) {

    seqs <- Biostrings::readDNAStringSet(seqfile)
  
    # TODO: 'trainingSet' name is expected, cannot vary
    load(db)
    
    cat("\nClassifying", length(seqs), "sequences from", seqfile,
        "against", length(trainingSet$taxonomy), "reference sequences",
        "at confidence >=", threshold, "using", processors, "processors.\n",
        file=stderr())

    if (!is.null(rand_seed)) {
      cat("Setting random seed:", rand_seed, "\n", file=stderr())
      set.seed(rand_seed)
    }
    ids <- DECIPHER::IdTaxa(
        seqs,
        trainingSet,
        strand="both",
        ...
    )

    l = make_lineage(ids, total_ranks=max(trainingSet$levels))
    out = data.frame(
        `Feature ID`=names(ids),
        Taxon=names(l),
        Confidence=round(unname(l), 1),
        check.names=F,
        stringsAsFactors=F
    )
    write.table(out, gzfile(tax_out), sep="\t", na="", quote=F, row.names=F)    
}

log = file(snakemake@log[[1]])
sink(log, type="output")
sink(log, type="message")

params=snakemake@params$par

idtaxa_assign(
    snakemake@input$seq,
    snakemake@input$db, 
    snakemake@output$tax,
    processors=snakemake@threads,
    threshold=params$confidence, # 60 (cautious) or 50 (sensible)
    bootstraps=get_default(params, "bootstraps", 100),
    rand_seed=get_default(params, "rand_seed", NULL)
)
