

#' Reads the amplicon pipeline results from a directory
#' 
#' @param ...: Additional parameters are passed on to `read_qiime_taxonomy()`
read_pipeline_results = function(INPUT,
                                 OTUTAB = file.path(INPUT, 'denoised_otutab.txt.gz'),
                                 TAX_DIR = file.path(INPUT, 'taxonomy'),
                                 OTU_SEQS = file.path(INPUT, 'denoised.fasta'),
                                 CMP_DIR = file.path(INPUT, 'cmp'),
                                 TREE_FILE = file.path(INPUT, 'denoised_tree.tre'),
                                 ITSX_PREFIX = file.path(INPUT, 'ITSx', 'out'),
                                 verbose = F,
                                 ...) {
  out = list()
  
  maybe_gz = function(f, stop=NULL) {
    if (file.exists(f))
      return(f)
    f = gsub('\\.gz(ip)?$', '', f)
    if (file.exists(f))
      return(f)
    if (!is.null(stop))
      stop(stop)
  }
  
  # OTU table
  if (verbose)
    cat('Reading OTU table...\n', file = stderr())
  OTUTAB = maybe_gz(OTUTAB, stop='OTU table not found.')
  out$otutab = read_otutab(OTUTAB)
  
  # read taxonomy
  if (verbose)
    cat('Reading OTU taxonomy...\n', file = stderr())
  out$taxonomy = list()
  for (tax_file in list.files(TAX_DIR, '\\.txt(\\.gz)?')) {
    tax_file = maybe_gz(file.path(TAX_DIR, tax_file))
    stopifnot(!is.null(tax_file))
    name = basename(gsub('\\.txt(\\.gz)?$', '', tax_file))
    s = strsplit(name, '-')[[1]]
    if (length(s) == 3) {
      # remove middle part of name, which is there for technical reasons and not relevant
      name = paste(s[1], s[3], sep = '-')
    }
    if (verbose)
      cat(paste0('...', name, '\n'), file = stderr())
    if (!is.null(out$taxonomy[[name]])) {
      warning(sprintf('Taxonomy "%s is present both in gzip-compressed and uncompressed form, used "%s"', name, tax_file))
    }
    out$taxonomy[[name]] = read_qiime_taxonomy(tax_file, ...)
  }

  # OTU sequences
  if (verbose)
    cat('Reading OTU sequences...\n', file = stderr())
  OTU_SEQS = maybe_gz(OTU_SEQS, stop='Denoised sequences not found')
  out$refseq = Biostrings::readDNAStringSet(OTU_SEQS)
  
  # Phylogenetic tree
  TREE_FILE = maybe_gz(TREE_FILE)
  if (!is.null(TREE_FILE)) {
    if (verbose)
      cat('Reading phylogenetic tree...\n', file = stderr())
    out$tree = ape::ladderize(ape::read.tree(TREE_FILE))
  }
  
  # comparisons to OTU sequences
  cmp_files = list.files(CMP_DIR, '\\.txt(\\.gz)?$')
  out$cmp = list()
  if (length(cmp_files) > 0) {
    if (verbose)
      cat(paste('Comparisons to OTU sequences found:',
                paste(cmp_files, collapse=', '), '\n'), file = stderr())
    for (f in cmp_files) {
      f = maybe_gz(file.path(CMP_DIR, f))
      name = gsub('\\.txt$', '', basename(f))
      out$cmp[[name]] = read.delim(f, header=F,
                                   col.names=c('query', 'target', 'ident'),
                                   colClasses=c('character', 'character', 'numeric'))
    }
  }

  # ITSx
  if (dir.exists(dirname(ITSX_PREFIX))) {
    if (verbose)
      cat('ITSx directory found, reading...\n', file = stderr())
    posfile = paste0(ITSX_PREFIX, '.positions.txt.gz')
    posfile = maybe_gz(posfile, stop=paste('ITSx positions file not found: ', posfile))
    out$itsx_results = read_itsx_pos(posfile)
    nd = paste0(ITSX_PREFIX, '_no_detections.txt')
    if (!file.exists(nd)) {
      # if file is missing, then there should be no non-detected sequences
      stopifnot(file.info(paste0(ITSX_PREFIX, '_no_detections.fasta'))$size == 0)
    } else {
      no_detections = read.delim(nd, stringsAsFactors=F, header=F)[, 1]
      out$itsx_results = rbind(
        out$itsx_results,
        data.frame(OTU=no_detections, SSU=NA, ITS1=NA, ITS2=NA, LSU=NA, ITSx_cat='no_ITS', stringsAsFactors=F)
      )
    }
    otus_only = setdiff(out$itsx_results$OTU, rownames(out$otutab))
    itsx_only = setdiff(rownames(out$otutab), out$itsx_results$OTU)
    if (length(c(otus_only, itsx_only)) > 0) {
      stop(sprintf('Mismatch between taxa in ITSx directory, some files may not be up to date. Missing in ITSx: %s. Missing in OTU table: %s',
                   paste(otus_only, collapse=','), paste(itsx_only, collapse=',')))
    }
  }
  
  out
}


print_list = function(x, n) {
  ifelse(length(x) > n, paste0(paste(x[1:n], collapse = ', '), '…'), paste(x, collapse =
                                                                             ', '))
}

#' Combines different data sources into a phyloseq object.
#' In contrast to the `phyloseq::phyloseq` function, it also
#' checks whether the OTU/sample names and numbers match, and
#' reports problems as warnings or errors.
#' 
#' The input can be either all the listed parameters, or a named list 
#' containing them.
#'
#' @param otutab: count matrix with row names set to OTUs and column names to samples
#' @param taxonomy: data frame with an OTU column, or list of data frames [first in list will be chosen]
#' @param meta: data frame of metadata, with 'sample' column or row names set to samples
#' @param tree: (optional) tree object, see `ape::read.tree`
#' @param refseq: (optional) `Biostrings::XStringSet` of OTU sequences
#' @param itsx_results: (optional) data frame with at least an OTU and ITSx_cat column, 
#'     will be added to the right of the taxonomy table
#
# It is possible to specify a named list containing all these objects as first arguments.
# All further provided arguments, if named according to one above, will override the data in this list.
make_phyloseq = function(...) {
  library(phyloseq)
  data = list(...)
  
  get_data = function(data, name, required = T) {
    if (!(name %in% names(data))) {
      if (name %in% names(data[[1]])) {
        return(data[[1]][[name]])
      } else if (required) {
        stop(sprintf("'%s' not supplied to make_phyloseq", name))
      } else {
        return(NULL)
      }
    } else {
      return(data[[name]])
    }
  }
  
  # OTUs
  otutab = as.matrix(get_data(data, 'otutab'))
  
  # metadata
  meta = get_data(data, 'meta')
  meta = as.data.frame(meta)
  if (!is.null(meta$sample)) {
    if (any(is.na(meta$sample)))
      stop("Found 'sample' column in metadata, but it has missing entries.\n")
    if (any(duplicated(meta$sample)))
      stop(
        sprintf(
          "Found 'sample' column in metadata, but it has duplicate entries: %s\n",
          print_list(meta$sample[duplicated(meta$sample)], 10)
        )
      )
    row.names(meta) = meta$sample
  }
  
  # compare sample names
  otus_only = setdiff(colnames(otutab), rownames(meta))
  meta_only = setdiff(rownames(meta), colnames(otutab))
  if (length(otus_only) > 0) {
    warning(
      sprintf(
        'Some samples were only found in OTU table, but not the metadata: %s\n',
        paste(otus_only, collapse = ', ')
      )
    )
  }
  if (length(meta_only) > 0) {
    warning(
      sprintf(
        'Some samples were only found in metadata, but not the OTU table: %s\n',
        paste(meta_only, collapse = ', ')
      )
    )
  }
  
  # taxonomy
  tax = get_data(data, 'taxonomy')
  if ('list' %in% class(tax)) {
    warning(
      sprintf(
        'taxonomy supplied as list to make_phyloseq, using first item (%s)\n',
        names(tax)[1]
      )
    )
    tax = tax[[1]]
    stopifnot(class(tax) %in% c('data.frame'))
  }
  stopifnot(!is.null(tax$OTU))
  # obtain taxonomic ranks
  ranks = attr(tax, 'ranks')
  if (is.null(ranks)) {
    ranks = get_data(data, 'tax_ranks')
    if (is.null(ranks)) {
      stop("No 'ranks' attribute found in taxonomy. Please specify taxonomic ranks with a 'tax_ranks' argument.")
    }
  }
  
  # compare OTU tab and taxonomy
  otus_only = setdiff(rownames(otutab), tax$OTU)
  tax_only = setdiff(tax$OTU, rownames(otutab))
  if (length(otus_only) > 0) {
    warning(sprintf(
      'Some OTUs were only found in OTU table, but not the taxonomy: %s\n',
      paste(otus_only, collapse = ', ')
    ))
  }
  if (length(tax_only) > 0) {
    warning(sprintf(
      'Some OTUs were only found in taxonomy, but not the OTU table: %s\n',
      paste(tax_only, collapse = ', ')
    ))
  }
  
  # add optional ITSx data
  itsx = get_data(data, 'itsx_results', required = F)
  if (!is.null(itsx)) {
    n = nrow(tax)
    tax = merge(tax, itsx[,c('OTU', 'ITSx_cat')], by='OTU', all.x=T)
    stopifnot(nrow(tax) == n)
    if (any(is.na(tax$ITSx_cat))) {
      stop(
        sprintf(
          'ITSx category not defined for some OTUs, is there an OTU label mismatch? OTUs: %s',
          paste(tax$OTU[is.na(tax$ITSx_cat)], collapse = ', ')
        )
      )
    }
  }
  
  tax = as.matrix(tax)
  if ('OTU' %in% colnames(tax)) {
    rownames(tax) = tax[, 'OTU']
    tax = tax[, colnames(tax) != 'OTU']
  }
  # make sure that the taxonomic ranks come first, otherwise, tax_glom will fail
  col_order = c(intersect(ranks, colnames(tax)), setdiff(colnames(tax), ranks))
  tax = tax[, col_order]
  
  physeq = phyloseq::phyloseq(otu_table(otutab, taxa_are_rows = T),
                              tax_table(tax),
                              sample_data(meta))
  
  # add optional data
  tree = get_data(data, 'tree', required = F)
  if (!is.null(tree)) {
    physeq2 = physeq
    phy_tree(physeq2) = tree
    if (phyloseq::ntaxa(physeq2) < phyloseq::ntaxa(physeq)) {
      warning('The phylogenetic tree contains less OTUs than the OTU table. Not added.\n')
    } else {
      physeq = physeq2
    }
  }
  
  seq = get_data(data, 'refseq', required = F)
  if (!is.null(seq)) {
    physeq2 = merge_phyloseq(physeq, refseq(seq))
    if (ntaxa(physeq2) < ntaxa(physeq)) {
      warning('There are less OTU sequences than OTUs in the count table. Not added.\n')
    } else {
      physeq = physeq2
    }
  }
  
  physeq
}


#' Reads QIIME-style OTU table
read_otutab = function(otu_file) {
  otutab = as.data.frame(suppressMessages(readr::read_tsv(otu_file)))
  row.names(otutab) = otutab[, 1]
  otutab = otutab[2:ncol(otutab)]
  as.matrix(otutab)
}


#' Reads a positions file from ITSx
read_itsx_pos = function(positions) {
  cols = c('OTU', 'len', 'SSU', 'ITS1', 'S58', 'ITS2', 'LSU', 'comment')
  domains = c('SSU', 'ITS1', 'S58', 'ITS2', 'LSU')
  x = suppressMessages(read.delim(
    positions,
    header=F,
    col.names=cols,
    stringsAsFactors=F
  ))
  x = cbind(
    x['OTU'], 
    apply(x[domains], 2, function(x) {
      out = sapply(strsplit(x, ': ', fixed=T), '[', 2)
      ifelse(out == 'Not found', NA, out)
    }),
    x['comment']
  )
  # remove --END--
  x = x[x$OTU != '--END--', ]
  x$ITSx_cat = with(x, ifelse(
    is.na(comment),
    'ok',
    gsub(' ', '_', gsub('! *$', '', gsub('only partial ', 'partial_', gsub('Broken or partial sequence, ', '', comment))))
  ))
  x[, c('OTU', 'SSU', 'ITS1', 'ITS2', 'LSU', 'ITSx_cat')]
}



#### Taxonomy ####

#' Reads taxonomic annotations in the tabbed format created by QIIME2
#' OTU IDs are in the first column
#' Missing values are automatically completed using `replace_missing_taxa`. 
#' Any additional arguments are forwarded to this function.
read_qiime_taxonomy = function(tax_file, ...) {
  tax = read.delim(tax_file,
                   colClasses='character',
                   na.strings=c('NA', 'Unassigned'),
                   stringsAsFactors=F)
  t = as.data.frame(expand_lineage(tax[, 2]), stringsAsFactors=F)
  ranks = names(t)
  t = cbind(OTU = tax[, 1], t)
  nranks = length(ranks)
  t = replace_missing_taxa(t, ranks, ...)
  t = as.data.frame(t, stringsAsFactors = F)
  attr(t, 'ranks') = ranks
  ## this code writes back the lineages, allowing to compare raw input with adjusted lineages
  # write.table(data.frame(
  #   `Feature ID`=t[, 'OTU'],
  #   Taxon=apply(sapply(ranks, function(r) paste0(substr(r, 1, 1), '__', gsub(' ', '_', t[,r]))), 1, function(x) paste(x, collapse='; ')),
  #   Confidence=tax[,3]
  # ), paste0(tax_file, '.write_back.txt'), sep='\t', quote=F, row.names=F)
  t
}


#' Reads assignment results from PROTAX-Fungi (https://github.com/psomervuo/protaxfungi).
#' Tested with runs downloaded from the PlutoF workbench (https://plutof.ut.ee/,
#' Laboratories -> Analysis Lab -> Sequence analysis -> New -> PROTAX-Fungi).
read_protax = function(protax_dir,
                       ranks = c('kingdom',
                                 'phylum',
                                 'class',
                                 'order',
                                 'family',
                                 'genus',
                                 'species'),
                       threshold = 0.9,
                       ...)
{
  # read probabilities
  protax_l = lapply(2:7, function(i) {
    f = file.path(protax_dir, sprintf('query%d.nameprob', i))
    lines = strsplit(readChar(f, file.info(f)$size), '\n', fixed = T)[[1]]
    fields = strsplit(lines, '\t', fixed = T)
    as.data.frame(t(sapply(fields, function(x)
      x[1:3])))
  })
  # OTU order should be same in all files
  for (d in protax_l) {
    stopifnot(d$OTU == protax_l[[1]]$OTU)
  }
  # matrix of probabilities
  probs = sapply(protax_l, function(d)
    as.numeric(d[[3]]))
  selected = probs >= threshold
  lineages = sapply(protax_l, function(d)
    d[[2]])
  # assemble lineages
  protax = t(sapply(1:nrow(probs), function(i) {
    pr = selected[i, ]
    if (!pr[1] | is.na(pr[1])) {
      return(rep(NA, 7))
    }
    # up to which level are probabilities equal or greater than threshold?
    lin = strsplit(lineages[i, max(which(pr))], ',', fixed = T)[[1]]
    c(lin, rep(NA, 7 - length(lin)))
  }))
  colnames(protax) = ranks
  rownames(protax) = protax_l[[1]][[1]]
  protax[protax == 'unk'] = NA
  # remove number at end of genus / species names
  protax = gsub('_[0-9]+$', '', protax, perl = T)
  out = replace_missing_taxa(protax, ranks = ranks, ...)
  out = cbind(data.frame(OTU = row.names(out)), as(out, 'data.frame'))
  attr(out, 'ranks') = ranks
  out
}

#' Expand QIIME-like taxonomic annotations into a matrix, replacing the prefixes
#' with known full rank names
expand_lineage = function(x,
                          pattern = ' *([a-z])__([^;]+)',
                          delim = ';',
                          # replace names starting with the same prefixes with the whole names taken from this vector
                          replace_ranks = c('kingdom',
                                            'phylum',
                                            'class',
                                            'order',
                                            'family',
                                            'genus',
                                            'species')) {
  ranks = strsplit(gsub(pattern, '\\1', x, perl = T), delim)
  lineages = strsplit(gsub(pattern, '\\2', x, perl = T), delim)
  named = lapply(1:length(ranks), function(i) {
    l = lineages[[i]]
    names(l) = ranks[[i]]
    l
  })
  rankcols = na.omit(unique(unlist(ranks)))
  cols = lapply(rankcols, function(r)
    sapply(named, '[', r))
  t = do.call(cbind, cols)
  # try replacing names
  i = sapply(rankcols, function(col)
    which(startsWith(replace_ranks, col))[1])
  rankcols[!is.na(i)] = replace_ranks[!is.na(i)]
  colnames(t) = rankcols
  rownames(t) = names(x)
  t
}


#' Replaces NA values in a taxonomy data frame with appropriate names based
#' on the names of the lower rank. In addition, adds a `def_rank` column
#' with the rank, at which an assignment is defined.
#' 
#' By default, NAs in species are converted to 'Genus sp.' or (e.g.) to
#' 'Phylum clone'. But the defaults can be changed (`unknown_species_fmt`).
#' NA in other ranks are converted to <lower_rank_name>_<rank>_unknown, e.g.
#' Eurotiomycetes_gen_unknown or Viridiplantae_ord_unknown. The behaviour
#' can be changed with `unknown_fmt`
#'
#' @param tax_tab: Taxonomy table (matrix or data frame)
#' @param ranks: names/indices of columns to consider (other will be left untouched)
#'    default: all columns
#' @param top_unknown: use the given name if the leftmost column (top rank) is `NA`
#' @param unknown_species_fmt: name to use for unknown species, when the
#'    genus is *not* known. This *must* be a format string (see `unknown_fmt`),
#'    where '%s' is replaced by the highest known rank name.
#' @param unknown_fmt: for all other ranks except for toplevel and species,
#'    the name is formed based on the highest known rank. This *must* be a format
#'    string with two replacement variables '%s'. The first is replaced by the 
#'    highest known rank name, the second by the first three letters of the rank.
#'    default: '%s_%s_unknown' (yielding e.g. Eurotiomycetes_gen_unknown).
#     To exchange their order, use something like '%2$s_%1$s' (gives 'gen_Eurotiomycetes').
#' @param species_name: name of species column.
#'   If present, underscores are replaced with spaces in species names,
#'   and NA in species column will be converted to 'Genus sp.' or (e.g.) to
#'   'Phylum clone', at least by default (see also `unknown_species_fmt`)
#' @param genus_name: If the genus name is not 'genus' and 'species_name'
#'   is present, this must be specified. Needed for formatting species names
#' @param not_def_rank: Value to fill in `def_rank` if the whole lineage 
#'   consists of `NA`s
replace_missing_taxa = function(tax_tab,
                           ranks = colnames(tax_tab),
                           top_unknown = 'Unknown',
                           unknown_species_fmt = '%s clone',
                           unknown_fmt = '%s_%s_unknown',
                           genus_name = ifelse('genus' %in% ranks, 'genus', NULL),
                           species_name = ifelse('species' %in% ranks, 'species', NULL),
                           not_def_rank = '<none>') {
  stopifnot(!is.null(top_unknown))
  # prepare species completion
  do_fill_spec = F
  if (!is.null(species_name) & !is.null(genus_name)) {
    stopifnot(species_name %in% ranks)
    stopifnot(genus_name %in% ranks)
    if (length(regmatches(unknown_species_fmt, gregexpr('%s', unknown_species_fmt, fixed=T))[[1]]) != 1)
      stop('replace_missing_taxa: unknown_species_fmt needs exactly one "%s" wildcard')
    # normalize species names: replace underscores by spaces
    tax_tab[, species_name] = gsub('_', ' ', tax_tab[, species_name], fixed=T)
    fill_ranks = ranks[ranks != species_name[1]]
    do_fill_spec = T
  } else {
    fill_ranks = ranks
  }
  # prepare filling of other ranks
  if (length(regmatches(unknown_fmt, gregexpr('%(\\d+\\$)?s', unknown_fmt, perl=T))[[1]]) != 2)
    stop('replace_missing_taxa: unknown_fmt needs exactly two "%s" wildcards')
  short_fill_ranks = substr(fill_ranks, 1, 3)
  n_fill = length(fill_ranks)
  # iterate through taxonomy
  # note: this may not be the fastest running code, but it was rather 
  # programmed to be understandable
  tax = t(apply(tax_tab, 1, function(x) {
    x_sel = x[ranks]
    defined = !is.na(x_sel)
    # fill toplevel rank if not defined
    if (!defined[1]) {
      x[ranks[1]] = x_sel[1] = top_unknown
      defined[1] = T
      def_i = 1
      def_rank = not_def_rank
    } else {
      def_i = max(which(defined))
      def_rank = ranks[def_i]
    }
    # fill species name
    if (do_fill_spec && !defined[species_name]) {
      fmt = if (defined[genus_name]) '%s sp.' else unknown_species_fmt
      x[species_name] = x_sel[species_name] = sprintf(fmt, x_sel[def_i])
    }
    # fill other ranks
    fill_lineage = x_sel[fill_ranks]
    defined = !is.na(fill_lineage)
    stopifnot(defined[1])  # toplevel rank should always be filled
    def_i = max(which(defined))
    if (any(!defined)) {
      replacements = sprintf(unknown_fmt, fill_lineage[def_i], short_fill_ranks)
      x[fill_ranks][!defined] = replacements[!defined]
    }
    c(x, def_rank = def_rank)
  }))
  tax = as.data.frame(tax, stringsAsFactors = F)
  tax$def_rank = factor(tax$def_rank, c(ranks, not_def_rank))
  tax
}

