library(stringr)
library(yaml) 


get_default <- function(list, item, default=NULL) {
    value <- list[[item]]
    if (is.null(value)) default else value
}

make_lineages <- function(headers) {
  desc <- str_split_fixed(headers, " ", 2)[,2]
  nranks <- max(str_count(desc, ";")) + 1
  tax <- str_split_fixed(desc, ";", nranks)
  stopifnot(!any(is.na(tax)))  # no rank should be undefined in input
  # separate ranks from names and ensure that the ranks are the same across
  # columns
  rank_names <- apply(tax, 2, function(x) str_split_fixed(x, "__", 2), simplify=F)
  stopifnot(sapply(rank_names, function(x) length(unique(x[,1]))) == 1)
  ranks <- sapply(rank_names, '[', 1)
  taxa_out <- sapply(rank_names, function(x) x[,2])
  # replace empty ranks with NAs: since internal undefined ranks should not be
  # empty (rank propagated earlier in pipeline),
  # we should not obtain any NAs in output lineages
  taxa_out[taxa_out == ""] <- NA
  # add root
  taxa_out <- cbind("Root", taxa_out)
  ranks <- c("rootrank", ranks)
  # get defined rank
  def_rank <- apply(taxa_out, 1, function(x) max(which(!is.na(x))))
  # assemble lineages
  lineages <- sapply(seq_along(def_rank), function(i) {
    d <- def_rank[i]
    x <- taxa_out[i,]
    paste(x[1:d], collapse="; ")
  })
  # ensure that there were no NAs
  # TODO: there can in theory be false positives here
  stopifnot(!grepl("; NA(;|$)", lineages, perl=T))
  # return all the info
  list(
    ranks = ranks,
    def_rank = def_rank,
    taxa = lineages
  )
}


idtaxa_train <- function(taxdb_file, trained_out,
                          processors=1,
                          max_group_size=NULL,
                          max_iterations=3,
                          allow_group_removal=F,
                          rand_seed=NULL) {
  # read data and make sure it's all in the same orientation
  cat("Reading input sequences from", "taxdb_file", "...\n", file=stderr())
  seqs <- Biostrings::readDNAStringSet(taxdb_file)
  
  # subsample: only for testing!!
  # seqs <- seqs[sample(length(seqs), 2000)]
  
  cat("Running 'OrientNucleotides'...\n", file=stderr())
  seqs <- DECIPHER::OrientNucleotides(seqs, processors=processors)

  # obtain lineages formatted as expected by IDTAXA
  seq_headers <- names(seqs)
  tax <- make_lineages(seq_headers)
  names(seqs) <- tax$taxa
  
  # vector for removed sequences (used at several places)
  remove <- logical(length(seqs))
  # taxa groups
  groups <- names(seqs)
  group_counts <- table(groups)
  u_groups <- names(group_counts)
  
  if (!is.null(rand_seed)) {
    cat("Setting random seed for pruning + LearnTaxa:", rand_seed, "\n", file=stderr())
    set.seed(rand_seed)
  }
  
  # prune the training set
  if (!is.null(max_group_size)) {
    for (i in which(group_counts > max_group_size)) {
      index <- which(groups == u_groups[i])
      keep <- sample(length(index), max_group_size)
      remove[index[-keep]] <- TRUE
    }
    rm_freqs <- sort(table(groups[remove]), decreasing=T)
    cat(sprintf(
      "Removed %d of %d sequences from large groups with > %d sequences:\n%s\n\n",
      sum(remove), length(groups), max_group_size,
      paste(paste(names(rm_freqs), rm_freqs, sep=': '), collapse='\n')
    ), file=stderr())
  }
  
  # iterative training
  # (this code is taken unmodified from the tutorial
  # http://www2.decipher.codes/Documentation/Documentation-ClassifySequences.html)
  probSeqsPrev <- integer() # suspected problem sequences from prior iteration
  for (i in seq_len(max_iterations)) {
    cat("Training iteration: ", i, "\n", sep="")
    # train the classifier
    trainingSet <- DECIPHER::LearnTaxa(
        seqs[!remove],
        names(seqs)[!remove]
    )
    # look for problem sequences
    probSeqs <- trainingSet$problemSequences$Index
    if (length(probSeqs) == 0) {
      cat("No problem sequences remaining.\n")
      break
    } else if (length(probSeqs) == length(probSeqsPrev) &&
               all(probSeqsPrev==probSeqs)) {
      cat("Iterations converged.\n")
      break
    }
    cat("Problematic sequences in iteration", i, ":\n",
        paste(seq_headers[!remove][probSeqs], collapse="\n"),
        "\n\n", file=stderr())
    
    if (i == max_iterations)
      break
    probSeqsPrev <- probSeqs
    # remove any problem sequences
    index <- which(!remove)[probSeqs]
    remove[index] <- TRUE  # remove all problem sequences
    if (!allow_group_removal) {
      # replace any removed groups
      missing <- !(u_groups %in% groups[!remove])
      missing <- u_groups[missing]
      if (length(missing) > 0) {
        index <- index[groups[index] %in% missing]
        remove[index] <- FALSE # don't remove
      }
    }
  }

  # # info plot
  # pdf(info_out, width=9, height=7.4)
  # plot(trainingSet)
  # dev.off()
  
  # we did not supply the 'rank' argument, instead we just set the ranks manually here
  trainingSet$ranks <- tax$ranks[trainingSet$levels]
  
  # finally save the databaes
  save(trainingSet, file=trained_out)
}

log = file(snakemake@log[[1]])
sink(log, type="output")
sink(log, type="message")

params <- read_yaml(snakemake@input$params)

idtaxa_train(
    snakemake@input$seq,
    snakemake@output$db, 
    #processors=snakemake@threads,
    max_group_size=get_default(params, "max_group_size", NULL),
    max_iterations=get_default(params, "max_iterations", 3),
    allow_group_removal=get_default(params, "allow_group_removal", F),
    rand_seed=get_default(params, "rand_seed", NULL)
)
