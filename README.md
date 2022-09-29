## Cancer Testis Antigen (CTA) mining using public gene expression data across 33 cancers
In this repo, we carry out a large-scale systematic analysis of the expressions of approx. 250 cancer testis antigens (CTAs) across 33 cancers, with the goal of mining for the top CTA candidates that can be used for new cancer therapies (targets for drug design).  The analysis also takes into account a standard measure of mutational load associated with each patient sample, known as the tumor mutation burden (TMB). This measure is known to positively correlate with the expression of certain CTAs. Therefore, this aspect is also looked at through this analysis to mine the best CTA candidates that are highly expressed in many cancers but have little to no expression in normal tissues.

### What are CTAs?
Cancer testis antigens (CTA) are a group of genes (with the end products being proteins) expressed in many human tumors but not in normal tissues except for testis and placenta. This tumor-specific pattern of expression, along with strong immunogenicity, identifies CTAs as ideal targets for tumor-specific immunotherapeutic approaches. 

CTA expressions have been found in melanoma, liver cancer, lung cancer, bladder cancer, and pediatric tumors such as neuroblastoma. Several clinical trials of CTA‐based vaccine therapy have been developed or are being currently run. So far, at least 140 CTAs have been identified, most of which are expressed during spermatogenesis. 

The Ludwig Institute for Cancer Research (LICR) maintains the ["CTDatabase",](http://www.cta.lncc.br/) which is a comprehensive reference of known CTAs and is continually being updated. For this analysis, the list of CTAs was obtained from this database. 

For more info on CTAs, see here:
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5528287/
- https://en.wikipedia.org/wiki/Cancer/testis_antigens

### Why exactly do we insist on finding CTAs that are expressed in cancers but not in normal tissues? In layman terms, why is this important for effective immunotherapy?
Due to the COVID pandemic, there is a lot more awareness and understanding of how vaccines work in general. In simple terms, immunotherapy can be understood as trying to make vaccines for cancer treatment or, in other words, trying to stimulate a cancer patient's own immune system to search for, attack, and kill cancer cells. Therefore, we need to identify antigens that are selectively expressed (or present) only in tumors and are not expressed much in normal tissues. If we are able to identify such antigens, then we can target those antigens to design vaccines that can then be used to train our own immune systems to identify and kill cancer cells very precisely, without damaging normal and healthy cells in our body. This is the basic idea, but of course, like anything else in biology, the actual mechanisms are much more complicated, and hence the need for rigorous research in the field.

### The data
All data used in the analysis are publicly available. The main gene expression matrix for the cancers along with relevant clinical data was obtained from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga), which is a very useful large-scale consortium effort that has a rich set of data across 33 cancers gathered systematically from 11000 patient samples. 

To compare the cancer gene expressions with their expressions in normal tissues, the normal human organ and tissues data from [The Genotype-Tissue Expression (GTEx) project.](https://gtexportal.org/home/) was downloaded.

Instead of downloading separately from the individual sources, a combined and convenient matrix format containing both the TCGA cancer and the GTEx normal expression data is available from the [UCSC Xena Browser.](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) This source also has the additional advantage that the cancer and the normal expression matrices have been reprocessed to normalize for batch effects. In other words, both the cancer and normal expressions have been integrated in an organized fashion, ready to be analyzed together.

The above matrix from Xena is about 7.8 GB in size and has 58581 rows (genes) and 19121 columns (patient samples). The matrix also has samples from a pediatric cancer dataset called [TARGET](https://ocg.cancer.gov/programs/target). This TARGET set of pediatric patients is out of the scope of this analysis, and hence, we exclude those columns. For our analysis here, we retain only the GTEx normal samples and TCGA adult cancer samples (18316 columns in total).


