#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse)
library(openxlsx)

# read mock frequencies
freq = read_tsv('mock/mock_data.txt') %>% 
  select(species=name, starts_with('rel_abund_')) %>% 
  pivot_longer(-species, names_to='sample', values_to='rel_mixed') %>% 
  mutate(sample = gsub('rel_abund_m', 'mock', sample))

dirs = list.dirs('results')
dirs = dirs[grepl('^results/.+?/workflow_.+?/.+?/[^/]+?__[^/]+?\\.{3,3}[^/]+?$', dirs, perl=T)]

if (length(dirs) > 0) {
  primer_combs = basename(dirs)
  s = split(dirs, primer_combs)

  for (primer_comb in names(s)) {
    outdir = file.path('mock_cmp', primer_comb)
    primer_dirs = s[[primer_comb]]
  # remove validation workflows, which don't have a sequence comparison file
    primer_dirs = primer_dirs[!grepl('_simple', primer_dirs, fixed=T)]
    seqs = file.path(primer_dirs, 'clusters.fasta')
    tabs = file.path(primer_dirs, 'otutab.txt.gz')
    mapping = file.path(primer_dirs, 'cmp', 'mock.txt')
    sel = file.exists(seqs)
    primer_dirs = primer_dirs[sel]
    tabs = tabs[sel]
    seqs = seqs[sel]
    stopifnot(file.exists(tabs) & file.exists(mapping))
    data = map_dfr(seq_along(primer_dirs), function(i) {
      # OTU tabs
      tab_file = tabs[i]
      groups = stringr::str_match(
        primer_dirs[i], 
        ".+?/(?<workflow>[^/]+)/workflow_.+?/(?<run>.+?)_(?<layout>[^/_]+)/(?<primers>[^/]+)"
      )[, 2:5, drop=F]
      d = suppressMessages(read_tsv(tab_file))
      names(d)[1] = 'OTU'
      # OTU -> mock reference mapping
      cmp = suppressMessages(read_tsv(
        mapping[i], 
        col_names=c('OTU', 'species', 'ident'),
        col_types='ccn'
      ))
      d = d %>% 
        pivot_longer(-OTU, names_to='sample', values_to='count') %>% 
        left_join(cmp, 'OTU') %>% 
        left_join(freq, c('species', 'sample'))
      bind_cols(groups, d)
    })

    data = data %>% 
      mutate(group=paste(workflow, run, layout)) %>% 
      group_by(group, sample, species, rel_mixed) %>% 
      summarise(
        count = sum(count), 
        n=n(), 
        OTU = paste(OTU, collapse=','), 
        .groups='drop'
      )
    
    leg = guide_legend(ncol=1)
    ggplot(data, aes(rel_mixed, count, colour=group, shape=group)) +
      facet_grid( ~ sample, scales='free') + 
      geom_point(aes(size=as.integer(as.factor(n))), pch=21, alpha=0.7) +
      geom_smooth(method='lm', formula='y~x', se=F, linewidth=0.3) +
      scale_size_continuous(range=c(1, 2)) +
      scale_colour_brewer(palette='Paired') +
      scale_x_log10() +
      scale_y_log10() +
      labs(x='Mixed relative concentration',
           y='Read count',
           colour='Workflow/run',
           shape='Workflow/run',
           size='Number of ASVs') +
      guides(colour=leg, shape=leg, size=leg) +
      theme_bw() +
      theme(strip.text.y=element_text(angle=0),
            legend.position='bottom')
    ggsave(file.path(outdir, 'mock.png'),
           width=16, height=12, units='cm', dpi=300)
    last_plot() + facet_grid(group ~ sample, scales='free')
    ggsave(file.path(outdir, 'mock_detail.png'),
                     width=16, height=25, units='cm', dpi=300)

    # create Excel file count tables
    wb = createWorkbook()
    hs = createStyle(border='Bottom')
    for (smp in split(data, data$sample)) {
      sample = smp$sample[1]
      m = smp %>% 
        select(-OTU, -sample, -n) %>% 
        mutate(species=if_else(is.na(species), "(other)", species)) %>% 
        arrange(-rel_mixed, group) %>% 
        pivot_wider(names_from=group, values_from=count, values_fill=0, names_sep=' ') %>% 
        column_to_rownames('species')
      # add totals on top
      m = rbind(total=colSums(m), m)
      # write data
      addWorksheet(wb, sample)
      writeData(wb, sample, c('sample', rownames(m)), startCol=1, headerStyle=hs)  # cluster names
      writeData(wb, sample, m, startCol=2, headerStyle=hs)  # counts
      freezePane(wb, sample, firstRow=T, firstCol=T)
      setColWidths(wb, sample, 1, 15)
      addStyle(wb, sample, hs, 2, 1:(ncol(m)+2))
      q = c(0, 0.5, 1)
      colors = c('#abdda4', '#ffffbf', '#fdae61')
      conditionalFormatting(wb, sample, cols=2, rows=3:(nrow(m)+1),
                            type='colourScale', style=colors,
                            rule=quantile(unlist(m[2:nrow(m), 1]), q, na.rm=T))
      conditionalFormatting(wb, sample, cols=3:(ncol(m)+2), rows=3:(nrow(m)+1),
                            type='colourScale', style=colors,
                            rule=quantile(unlist(m[2:nrow(m), 3:ncol(m)]), q, na.rm=T))
    }
    saveWorkbook(wb, file=file.path(outdir, 'counts.xlsx'), overwrite=T)
  }    
}
