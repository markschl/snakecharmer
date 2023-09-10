#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(dtplyr)
library(tidyverse, warn.conflicts=F)
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
               col_names=c('path', 'workflow', 'run', 'layout', 'primers', 'OTU', 'cluster', 'id'))
# similar clusters
self_cmp = read_tsv(file.path(pipeline_dir, 'cmp', 'self_cmp.txt'),
               col_names=c('primers', 'cluster', 'other', 'id'),
               col_types='cccd') %>% 
  lazy_dt() %>% 
  group_by(primers, cluster) %>% 
  arrange(-id) %>% 
  slice(1:min(n(), 5)) %>% 
  summarise(
    similar_clusters=paste(sprintf('%s (%.1f%%)', other, id), collapse=', '),
    .groups='drop'
  ) %>% 
  as_tibble() %>% 
  split(.$primers) %>% 
  map(~ select(.x, cluster, similar_clusters) %>% deframe())


s = split(cmp, cmp$primers)

for (primers in names(s)) {
  cmp_primer = s[[primers]]
  
  data = map_dfr(split(cmp_primer, cmp_primer$path), function(d) {
    tab_file = file.path(d$path[1], 'otutab.txt.gz')
    otus = suppressMessages(read_tsv(tab_file))
    names(otus)[1] = 'OTU'
    lazy_dt(otus) %>% 
      pivot_longer(-OTU, names_to='sample', values_to='count') %>% 
      left_join(d, 'OTU') %>% 
      group_by(workflow, run, layout, sample, cluster) %>% 
      summarise(count=sum(count), n=n()) %>% 
      as_tibble()
  })
  
  cluster_file = file.path(pipeline_dir, 'cmp', paste0('clusters.', primers, '.fasta'))
  seqs = Biostrings::readDNAStringSet(cluster_file)
  seqlen = tibble(cluster=names(seqs), seq_length=lengths(seqs))
  
  data = data %>% 
    left_join(seqlen, 'cluster')
  
  count_tabs = data %>% 
    split(.$sample) %>% 
    map(function(smp) {
      m = smp %>% 
        select(-n, -sample) %>% 
        arrange(workflow) %>% 
        pivot_wider(names_from=c(workflow, run, layout), values_from=count,
                    values_fill=0, names_sep=' ') %>% 
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
    m_sel = m[,2:ncol(m), drop=F]
    m_rel = as.data.frame(t(t(m_sel) / colSums(m_sel) * 2))
    for (c in names(m_rel)) class(m_rel[[c]]) = 'percentage'
    # write data
    addWorksheet(wb, sample)
    writeData(wb, sample, c('sample', rownames(m)), startCol=1, headerStyle=hs)  # cluster names
    writeData(wb, sample, c('similar', self_cmp[[primers]][rownames(m)]), startCol=2,
              headerStyle=hs)  # similar clusters
    writeData(wb, sample, m, startCol=3, headerStyle=hs)  # counts
    writeData(wb, sample, m_rel, startCol=ncol(m)+4, headerStyle=hs)  # relative counts
    freezePane(wb, sample, firstRow=T, firstCol=T)
    setColWidths(wb, sample, 2, 25)
    addStyle(wb, sample, hs, 2, 1:(3*ncol(m)+3))
    colors = c('#abdda4', '#ffffbf', '#fdae61')
    q = c(0, 0.5, 1)
    conditionalFormatting(wb, sample, cols=4:(ncol(m)+2), rows=3:(nrow(m)+1),
                          style=colors, rule=quantile(unlist(m[2:nrow(m),]), q, na.rm=T), 
                          type='colourScale')
    conditionalFormatting(wb, sample, cols=(ncol(m)+4):(2*ncol(m)+3), rows=3:(nrow(m)+1),
                          style=colors, rule=quantile(unlist(m_rel[2:nrow(m_rel),]), q, na.rm=T),
                          type='colourScale')
    colors = c('#ffb3b3', '#f7f7f7', '#80cdc1')
    q = c(0, 0.5, 1)
    conditionalFormatting(wb, sample, cols=3, rows=3:(nrow(m)+1),
                          style=colors, rule=quantile(m[2:nrow(m), 1], q, na.rm=T), 
                          type='colourScale')
  }
  saveWorkbook(wb, file=file.path(pipeline_dir, 'cmp', paste0(primers, '.xlsx')), 
               overwrite=T)
}

