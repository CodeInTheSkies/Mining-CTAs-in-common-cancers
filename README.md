## Cancer Testis Antigen (CTA) mining using public gene expression data across 33 cancers
In this repo, we carry out a large-scale systematic analysis of the expressions of approx. 250 cancer testis antigens (CTAs) across 33 cancers, with the goal of mining for the top CTA candidates that can be used for new cancer therapies (targets for drug design).  The analysis also takes into account a standard measure of mutational load associated with each patient sample, known as the tumor mutation burden (TMB). This measure is known to correlate positively with the expression of certain CTAs. Therefore, this aspect is also looked at through this analysis to mine the best possible CTAs.

### What are CTAs?
Cancer testis antigens (CTA) are a group of genes (with the end products being proteins) expressed in many human tumors but not in normal tissues except for testis and placenta. This tumor-specific pattern of expression, along with strong immunogenicity, identifies CTAs as ideal targets for tumor-specific immunotherapeutic approaches. 

CTA expressions have been found in melanoma, liver cancer, lung cancer, bladder cancer, and pediatric tumors such as neuroblastoma. Several clinical trials of CTA‐based vaccine therapy have been developed or are being currently run. So far, at least 140 CTAs have been identified, most of which are expressed during spermatogenesis. 

The Ludwig Institute for Cancer Research (LICR) maintains the ["CTDatabase",](http://www.cta.lncc.br/) which is a comprehensive reference of known CTAs and is continually being updated. For this analysis, the list of CTAs was obtained from this database. 

### The data
All data used in the analysis are publicly available. The main gene expression matrix for the cancers along with relevant clinical data was obtained from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga), which is a very useful large-scale consortium effort that has a rich set of data across 33 cancers gathered systematically from 11000 patient samples. 

To compare the cancer gene expressions with their expressions in normal tissues, we downloaded the normal tissue data for all human organs from [The Genotype-Tissue Expression (GTEx) project.](https://gtexportal.org/home/)

The most convenient and combined form of both the cancer and the normal expression data is available in a processed and ready-to-use format from the [UCSC Xena Browser.](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
