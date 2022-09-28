## Cancer Testis Antigen (CTA) mining using public gene expression data across 33 cancers
In this repo, we carry out a large-scale systematic analysis of the expressions of approx. 250 cancer testis antigens (CTAs) across 33 cancers, with the goal of mining for the top CTA candidates that can be used for new cancer therapies.  The analysis also takes into account a standard measure of mutational load associated with each patient sample, known as the tumor mutation burden (TMB).

### The data
All data used in the analysis are publicly available. The main gene expression matrix for the cancers along with relevant clinical data was obtained from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga), which is a very useful large-scale consortium effort that has a rich set of data across 33 cancers gathered systematically from 11000 patient samples. 

To compare the cancer gene expressions with their expressions in normal tissues, we downloaded the normal tissue data for all human organs from [The Genotype-Tissue Expression (GTEx) project.](https://gtexportal.org/home/)

The most convenient and combined form of both the cancer and the normal expression data is available in a processed and ready-to-use format from the [UCSC Xena Browser.](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
