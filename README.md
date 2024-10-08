# gro-crispri-ctcf

This repo contains the code to analyze and generate the plots seen in [Multi-locus CRISPRi targeting with a single truncated guide RNA](https://www.biorxiv.org/content/10.1101/2023.10.20.563306v1.article-metrics).
<br><br> To get started, clone the repo. There are jupyter notebooks that contain the code libraries and GEO data accessions to generate the plots.
Users will need to install jupyter notebook using `pip install notebook`. 
<br>

Each notebook contains install commands for R packages required for execution. 

## Clone repo
`gh repo clone broadinstitute/gro-crispri-ctcf`

## Descriptions
- `00_rnaseq_sg4_dream.ipynb` - RNA-seq analysis of sg4 in Jurkat cells
- `01_rnaseq_sg8_dream.ipynb` - RNA-seq analysis of sg8 in Jurkat cells
- `02_rnaseq_sgCD81_dream.ipynb` - RNA-seq analysis of sgCD81 in a375 cells
- `03-chipseq-sg4-ctcf-expected-sites.ipynb` - ChIP-seq CTCF analysis at sg4 perfect match sites in Jurkat cells
- `04-chipseq-sg4-h3k9me3-expected-sites.ipynb` - ChIP-seq H3K9me3 analysis at sg4 perfect match sites in Jurkat cells
- `05-chipseq-sg4-ctcf-dream.ipynb` - ChIP-seq CTCF differential analysis at ~880 putative CTCF binding sites in Jurkat cells using sg4
- `06-chipseq-sg4-ctcf-motif.ipynb` - ChIP-seq CTCF motif analysis at sg4 significant sites
- `07-chipseq-sg4-profile-maps.ipynb` - ChIP-seq CTCF, H3K9me3 and dCas9 profile heatmaps at sg4 perfect match sites
- `08-chipseq-sg8-ctcf-expected-sites.ipynb` - ChIP-seq CTCF analysis at sg8 perfect match sites in Jurkat cells
- `09-chipseq-sg8-h3k9me3-expected-sites.ipynb` - ChIP-seq H3K9me3 analysis at sg8 perfect match sites in Jurkat cells
- `10-chipseq-sg8-ctcf-dream.ipynb` - ChIP-seq CTCF differential analysis at ~880 putative CTCF binding sites in Jurkat cells using sg8
- `11-chipseq-sg8-ctcf-motif.ipynb` - ChIP-seq CTCF motif analysis at sg8 significant sites
- `12-chipseq-sg4-ctcf-partial-match-analysis.ipynb` - ChIP-seq CTCF motif analysis at sg4 partial match sites

## Accompanying scripts in R and python folders
- `filter_motif_by_pam.py` - Extract list of TF Motifs with NGG PAM
- `print_jaspar_motifs.R` - Create document of printed TF Motifs
- `tf-get-accession.py` - Extract ENCODE accessions given list of TFs
