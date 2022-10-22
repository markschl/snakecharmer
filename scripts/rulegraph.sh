#!/usr/bin/env bash

snakemake -d test --rulegraph \
    ITS taxonomy denoise config quality cmp \
    clean clean_all clean_taxdb clean_tax clean_cmp clean_itsx | 
    dot -Tpng > rulegraph.png
