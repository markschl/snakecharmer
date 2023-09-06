#!/usr/bin/env bash

snakemake -d test --rulegraph \
    config samples unique_samples quality \
    ITS taxonomy denoise cmp \
    clean clean_all clean_tax clean_cmp clean_itsx | 
    dot -Tpng > rulegraph.png
