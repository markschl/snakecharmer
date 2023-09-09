#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse)
library(data.table)
library(openxlsx)

args = commandArgs(trailingOnly=T)
if (length(args) == 0) {
  warning('No directory supplied, analyzing contents of results in currend directory')
  pipeline_dir = '.'
} else if (length(args)==1) {
  pipeline_dir = args[1]
}

# read results of sequence comparison
cmp = read_tsv(file.path(pipeline_dir, 'cmp', 'cmp.txt'),
               col_names=c('workflow', 'primers', 'layout', 'OTU', 'cluster', 'id'))
# similar clusters
self_cmp = read_tsv(file.path(pipeline_dir, 'cmp', 'self_cmp.txt'),
               col_names=c('primers', 'cluster', 'other', 'id'),
               col_types='cccd') %>% 
  group_by(primers, cluster) %>% 
  arrange(-id) %>% 
  summarise(similar_clusters=paste(sprintf('%s (%.1f%%)', other, id), collapse=', '), .groups='drop') %>% 
  split(.$primers) %>% 
  map(~ select(.x, cluster, similar_clusters) %>% deframe())

# get list of results directories
dirs = list.dirs(file.path(pipeline_dir, 'results'))
dirs = dirs[grepl('^results/.+?/workflow_.+?/.+?/[^/]+?__[^/]+?\\.{3,3}[^/]+?$', dirs, perl=T)]

if (length(dirs) > 0) {
  primer_combs = basename(dirname(dirs))
  s = split(dirs, primer_combs)
  
  for (primer_comb in names(s)) {
    primer_dirs = s[[primer_comb]]
    seqs = file.path(primer_dirs, 'clusters.fasta')
    tabs = file.path(primer_dirs, 'otutab.txt.gz')
    sel = file.exists(seqs)
    primer_dirs = primer_dirs[sel]
    tabs = tabs[sel]
    seqs = seqs[sel]
    stopifnot(file.exists(tabs))
    data = map_dfr(seq_along(primer_dirs), function(i) {
      # OTU tabs
      tab_file = tabs[i]
      workflow = basename(dirname(dirname(dirname(dirname(tab_file)))))
      layout = basename(dirname(tab_file))
      d = suppressMessages(read_tsv(tab_file))
      names(d)[1] = 'OTU'
      d %>% 
        pivot_longer(-OTU, names_to='sample', values_to='count') %>% 
        mutate(workflow=!!workflow,
               layout=!!layout)
    })
    
    stopifnot(!is.na(data$count))
    
    data = data %>% 
      # add comparison to clusters
      left_join(cmp, c('workflow', 'layout', 'OTU'))
    
    data = data.table(data)
    data = data[, .(count=sum(count), n=.N), by=.(workflow, layout, sample, cluster)]

    count_tabs = data %>% 
      split(.$sample) %>% 
      map(function(smp) {
        m = smp %>% 
          select(-n, -sample) %>% 
          arrange(workflow) %>% 
          pivot_wider(names_from=c(workflow, layout), values_from=count, values_fill=0, names_sep=' ') %>% 
          column_to_rownames('cluster')
      })
    
    # create Excel file with all tabs
    wb = createWorkbook()
    hs = createStyle(border='Bottom')
    for (sample in names(count_tabs)) {
      m = count_tabs[[sample]]
      m = m[rowSums(m) > 0,, drop=F]
      m = m[order(-rowSums(m)),,drop=F]
      # add totals on top
      m = rbind(total=colSums(m), m)
      # relative version
      m_rel = as.data.frame(t(t(m) / colSums(m) * 2))
      for (c in names(m_rel)) class(m_rel[[c]]) = 'percentage'
      # write data
      addWorksheet(wb, sample)
      writeData(wb, sample, c('sample', rownames(m)), startCol=1, headerStyle=hs)  # cluster names
      writeData(wb, sample, c('similar', self_cmp[[primer_comb]][rownames(m)]), startCol=2,
                headerStyle=hs)  # similar clusters
      writeData(wb, sample, m, startCol=3, headerStyle=hs)  # counts
      writeData(wb, sample, m_rel, startCol=ncol(m)+4, headerStyle=hs)  # relative counts
      freezePane(wb, sample, firstRow=T, firstCol=T)
      setColWidths(wb, sample, 2, 24)
      addStyle(wb, sample, hs, 2, 1:(3*ncol(m)+3))
      colors = c('#abdda4', '#ffffbf', '#fdae61')
      q = c(0, 0.5, 1)
      conditionalFormatting(wb, sample, cols=3:(ncol(m)+2), rows=3:(nrow(m)+1),
                            style=colors, rule=quantile(unlist(m[2:nrow(m),]), q, na.rm=T), 
                            type='colourScale')
      conditionalFormatting(wb, sample, cols=(ncol(m)+4):(2*ncol(m)+3), rows=3:(nrow(m)+1),
                            style=colors, rule=quantile(unlist(m_rel[2:nrow(m_rel),]), q, na.rm=T),
                            type='colourScale')
    }
    saveWorkbook(wb, file=file.path(pipeline_dir, 'cmp', paste0(primer_comb, '.xlsx')), 
                 overwrite=T)
  }
}
