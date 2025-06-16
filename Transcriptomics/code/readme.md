### File Names

naming scheme is `##-description` where \## is the number of the step (step 1 being ran first, step 2 second...) for code relevant to predominant workflow

### General Workflow

Retrieved egg transcripts STAR count matrix from Giles Goetz (NOAA).

-   `01-deseq2_with_STAR.Rmd` -- Pre-filters STAR count matrix and runs DESeq2 to determine differentially expressed genes (DEGs)

-   `02-blast_with_cds.Rmd` -- Runs BLAST using NCBI manila clam reference genome GCF_026571515.1 as query file.

-   `03-star_deg_annot.Rmd` -- Combines DEG results and BLAST results to get annotation information (SPID and GO terms). Also generates background (all expressed genes post filter) and gene (DEG) SPID lists for enrichment analysis via DAVID

-   `04-goslim.Rmd` -- Uses GO terms identified via BLAST to get 'slims' for all biological process GO terms.

-   `05-chisq_rnatypes.Rmd` -- Assesses if type distribution was different between control and primed eggs
