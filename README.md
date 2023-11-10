# gro-crispri-ctcf

This repo contains the code to analyze and generate the plots seen in "Multi-locus CRISPRi targeting with a single truncated guide RNA".
To get started, clone the repo. There are jupyter notebooks that contain the code libraries and GEO data accessions to generate the plots.

## clone repo
`gh repo clone broadinstitute/gro-crispri-ctcf`

## description
- 00_rnaseq_sg4_dream.ipynb - RNA-seq analysis of sg4 in Jurkat cells
- 01_rnaseq_sg8_dream.ipynb - RNA-seq analysis of sg8 in Jurkat cells
- 02_rnaseq_sgCD81_dream.ipynb - RNA-seq analysis of sgCD81 in a375 cells
- 03-chipseq-sg4-ctcf-expected-sites - ChIP-seq CTCF analysis at sg4 perfect match sites in Jurkat cells
- 04-chipseq-sg4-h3k9me3-expected-sites - ChIP-seq H3K9me3 analysis at sg4 perfect match sites in Jurkat cells
- 05-chipseq-sg4-ctcf-dream - ChIP-seq CTCF differential analysis at ~880 putative CTCF binding sites in Jurkat cells using sg4
- 06-chipseq-sg4-ctcf-motif - ChIP-seq CTCF motif analysis at sg4 significant sites
- 07-chipseq-sg4-profile-maps - ChIP-seq CTCF and H3K9me3 profile heatmaps at sg4 perfect match sites
- 08-chipseq-sg8-ctcf-expected-sites - ChIP-seq CTCF analysis at sg8 perfect match sites in Jurkat cells
- 09-chipseq-sg8-h3k9me3-expected-sites - ChIP-seq H3K9me3 analysis at sg8 perfect match sites in Jurkat cells
- 10-chipseq-sg8-ctcf-dream - ChIP-seq CTCF differential analysis at ~880 putative CTCF binding sites in Jurkat cells using sg8
- 11-chipseq-sg8-ctcf-motif - ChIP-seq CTCF motif analysis at sg8 significant sites
