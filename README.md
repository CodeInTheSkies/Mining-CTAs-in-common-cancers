## Cancer Testis Antigen (CTA) mining using public gene expression data across 33 cancers
In this repo, a large-scale systematic analysis of the expressions of approx. 250 cancer testis antigens (CTAs) across 33 cancers is carried out, with the goal of mining for the top CTA candidates that can be used for new cancer therapies (targets for drug design).  The analysis also takes into account a standard measure of mutational load associated with each patient sample, known as the tumor mutation burden (TMB). This measure is known to positively correlate with the expression of certain CTAs. Therefore, this aspect is also looked at through this analysis to mine the best CTA candidates that are highly expressed in many cancers but have little to no expression in normal tissues.

### What are CTAs?
Cancer testis antigens (CTAs) are a group of genes (with the end products being proteins) expressed in many human tumors but not in normal tissues except for testis and placenta. This tumor-specific pattern of expression, along with strong immunogenicity, identifies CTAs as ideal targets for tumor-specific immunotherapeutic approaches. 

CTA expressions have been found in melanoma, liver cancer, lung cancer, bladder cancer, and pediatric tumors such as neuroblastoma. Several clinical trials of CTA‐based vaccine therapy are being developed or are being currently run. So far, at least 140 CTAs have been identified, most of which are expressed during spermatogenesis. 

The Ludwig Institute for Cancer Research (LICR) maintains the ["CTDatabase",](http://www.cta.lncc.br/) which is a comprehensive reference of known CTAs and is continually being updated. For this analysis, the list of CTAs was obtained from this database. 

For more info on CTAs, see here:
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5528287/
- https://en.wikipedia.org/wiki/Cancer/testis_antigens

### Why exactly is finding CTAs expressed in cancers but not in normal tissues interesting? Why is this important for effective immunotherapy?
Due to the COVID pandemic, there is a lot more awareness and understanding of how vaccines work in general. In simple terms, immunotherapy can be understood as trying to make vaccines for cancer treatment or, in other words, trying to stimulate a cancer patient's own immune system to search for, attack, and kill cancer cells. Therefore, we need to identify antigens that are selectively expressed (or present) only in tumors and are not expressed much in normal tissues. If we are able to identify such antigens, then we can target those antigens to design vaccines that can then be used to train the immune system to identify and kill cancer cells very precisely, without damaging normal and healthy cells in the body. This is the basic idea, but of course, like anything else in biology, the actual mechanisms of how this could work are much more complicated, and hence the need for rigorous research in the field.

### The data
All data used in the analysis are publicly available. The main gene expression matrix for the cancers along with relevant clinical data was obtained from [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga), which is a very useful large-scale consortium effort that has a rich set of data across 33 cancers gathered systematically from 11000 patient samples. 

To compare the cancer gene expressions with their expressions in normal tissues, the normal human organ and tissue data from [The Genotype-Tissue Expression (GTEx) project](https://gtexportal.org/home/) was downloaded.

Instead of downloading separately from TCGA and GTEx, a combined and convenient matrix format containing both the TCGA cancer and the GTEx normal expression data is available from the [UCSC Xena Browser.](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) This source also has the additional advantage that the cancer and the normal expression matrices have been reprocessed to normalize for batch effects. In other words, both the cancer and normal expressions have been integrated in an organized fashion, ready to be analyzed together.

The above matrix from Xena is about 7.8 GB in size and has 58581 rows (genes) and 19121 columns (patient samples). The matrix also has samples from a pediatric cancer dataset called [TARGET](https://ocg.cancer.gov/programs/target). This TARGET set of pediatric patients is out of the scope of this analysis, and hence, we exclude those columns. For our analysis here, we retain only the GTEx normal samples and TCGA adult cancer samples (18316 columns in total).


The full list of CTAs chosen for analysis, a TCGA cancer name abbreviations lookup table, and some phenotype data needed for the analysis have been uploaded to the repo folder [smalldata](/smalldata).

For the TMB data, a single file containing TMB values for all cancer samples from TCGA could not be found for download. Instead, TMB values were downloadable as individual files for each cancer. These were downloaded, one per cancer, and then read one by one to populate a new TMB clinical data column.

All clinical data, including the TMB values for individual cancers, were downloaded from this [cBioportal page.](https://www.cbioportal.org/datasets)

From the above link, search for "pancan" ([The Pan-Cancer Atlas](https://gdc.cancer.gov/about-data/publications/pancanatlas)). All cancers available will be listed, and the links for each cancer will lead to its corresponding data.

For example, for lung cancer, the particular link would be:
[https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018](https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018)
From this main lung cancer page, choose the "clinical data" sub-tab and then choose "TMB(nonsynomymous)" and "Tumor Type" columns. All else can be deselected. Then download using the "down arrow" button below the "custom selection" button near the left of the bottom search bar (it should look like a cloud). This will download a tsv table as a separate, small-size text file with a few relevant clinical columns.

### The analysis
The analysis involves a number of data wrangling steps to clean the data and systematically match clinical data for each patient sample. Then the normal and cancer samples are binned into predefined groups including groups based on TMB ranges. For the final heatmaps (one per cancer), z-scores for the whole matrix need to be computed (just before plotting the heatmaps). Instead of the usual practice of row-wise z-scoring for gene-expression heatmaps, it was decided that computing z-scores for the entire matrix is a better method for this analysis as then the relative expression levels would have a global scale. This would be advantageous for bringing the most interesting CTAs to the top based on their expression values relative to the full set of CTAs in the matrix. Cancer and normal expressions can be fairly compared within a single heatmap if the z-scores are computed for the whole matrix.

The code here also supports subgrouping for certain cancers such as breast cancer, where it would be beneficial to divide the cancer samples into breast cancer subtypes. To achieve subgrouping, we use clinical data columns that have subtype information as needed depending on the type of cancer being analyzed.

The code is written such that we can analyze multiple cancers in a single run by specifying a list of abbreviations from the [TCGA abbreviations Table](/smalldata/tcga_abbr.txt) (see comments in the code for further details).

### Results
We now describe some sample results, interpretations, and follow-up analyses. 

Here is a sample heatmap section that has the top few and the bottom few CTAs for breast cancer, with patient samples grouped by subtypes and low, mid, and high TMB groups.

#### Top few and bottom few CTAs for Breast Cancer (with subtypes)
![Top few and bottom few CTAs for Breast Cancer (with subtypes)](/plots/HtMapTCGA_andNormals_CTAvsTMB_BreastCancer_Short.png)

A description of the elements of the plot is as follows. The heatmap has 4 main sections. The central section is the main heatmap with rows corresponding to the CTAs (genes) and columns corresponding to patient groups (normal tissues and cancer subgroups). The normal and the cancer portions are separated by a thin white vertical space. The bottom annotation is a column-wise bar plot with each bar representing the "cohort size" or the number of patients in each group. To the right, there are two special row annotations consisting of box and whisker plots. For each gene (or CTA), these boxplots represent that gene's average expression level (log2 scale, log2(normCount+1)) within the normal cohort excluding testis (left) and the cancer cohort (right), respectively. The center line in the boxes represents the median. The testis samples are excluded from the boxplot cohort because a high expression of CTAs in the testis is well known by definition, and so testis expression can be excluded to prevent unnecessary skewing of the median line towards high expression.

The vertical order of the CTAs is decided by the difference between the median value of the cancer boxplot and the median value of the normal boxplot. The heatmap is sorted in descending order of this difference from top to bottom. In other words, the heatmap is vertically sorted such that CTAs that have the maximum difference between the cancer median and the normal median would appear at the top.

#### Full-sized heatmap for Breast Cancer (all CTAs)
Here is the full-sized heatmap for Breast Cancer that has all CTAs analyzed. This is the full version of the above shortened heatmap.
![Fullsized plot](/plots/HtMapTCGA_andNormals_CTAvsTMB_BreastCancer_Full.png)

#### Full-sized heatmaps for many other cancers
A [large multipage PDF file](/plots/HtMapTCGA_andNormals_CTAvsTMB_Combined.pdf) containing full-sized heatmaps for a number of cancers analyzed thus far is available for download in [plots](/plots/HtMapTCGA_andNormals_CTAvsTMB_Combined.pdf).

### Follow-up analyses
