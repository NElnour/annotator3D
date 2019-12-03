# Three-dimensional Annotation of Motifs into Topologically-associating Domains

<!-- badges: start -->
<!-- badges: end -->

![](./inst/extdata/Elnour_N_A1.png)

Three-dimensional organization of the chromosome is a factor of gene expression regulation. Visualizations of 3D information is often presented separately from genome/epigenome features and motifs. One existing approach to integrate topological information is simply layering the TADs using GViz/ggbio. While powerful tools to correlate several annotation features, there is no current function to highlight the membership of features and motifs to distinct TADs using GViz/ggbio.

The goal of annotator3D is to annotate motifs into hierarchically clustered TADs. annotator3D does so by:
1. Soft-classifying motifs into TADs by positions (functions found in `classification.R`);
2. Extending GViz with a clusteredAnnotation class to highlight TAD membership (work in progress)

## Installation

You can install the released version of annotator3D with:

``` r
require(devtools)
devtools::install_github("NElnour/annotator3D", build_vignettes = TRUE)
library(annotator3D)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
data("mESChromInteractions")
data("matchedMotifs"")

tads <- cluster(mESChromInteractions, offset=40000)
```
![](./inst/extdata/dendrogram_tads.png)

``` r
tads.cl <- classify(matchedMotifs, tads, offset=40000)
head(tads.cl)

$H3K9me3_783_cluster_hc1_hc2_0.15_L2
 [1] "116: 0.02"     "114: 0.02"     "109-2-2: 0.04" "121: 0.02"     "122-2-1: 0.02" "131-1-2: 0.02" "137-2-2: 0.02"
 [8] "152-2-2: 0.06" "145: 0.02"     "144-1: 0.02"   "144-2: 0.04"   "140-2: 0.04"   "138-1: 0.02"   "141-1: 0.04"  
[15] "141-2-1: 0.02" "142-1: 0.02"   "143: 0.04"     "70-2: 0.02"    "65-1: 0.02"    "60: 0.02"      "50: 0.02"     
[22] "52: 0.02"      "56: 0.04"      "57: 0.02"      "46-1: 0.04"    "28-2: 0.02"    "22: 0.02"      "23: 0.02"     
[29] "17-1: 0.02"    "24: 0.02"      "10: 0.02"      "9: 0.02"       "12-2: 0.02"    "15: 0.04"      "16: 0.02"     
[36] "6-2-2-1: 0.02" "6-2-2-2: 0.02" "8: 0.02"       "3: 0.02"       "1-1: 0.02"    
```
Each number before the column corresponds to a TAD hierarchical cluster. The number after the colon corresponds to the cluster membership probability of the motif.

## Contributions
The author of this package is Nada Elnour. The functions within this package include:
1. `cluster`
    * `tad.group`
2. `classify`
    * `bin_motifs`
    * `calculate_overflow`
    
Wrapper function `cluster` relies on agglomerative clustering by `stats::hclust` function. It converts the resultant clusters into a dendrogram that are cut to generate predicted TADs and sub-TADs. It also uses tad.group and Bhel's `LinearizeNestedList` function to get the full path within sub-TADs from the root to the leaf. The user has the option to plot the predicted sub-TADs individually, and the output of `cluster` is:

1. a dataframe of hierarchically clustered chromsome segments, and their memeberships
2. a dendrogram plot whose colours denote the clusters.

Function `classify` relies on `bin_motifs` and `calculate_overflow` to calculate a soft classification score of motifs into the TADs and sub-TADs by position. It outputs a list of motifs the TADs/sub-TADs for which their scores are non-zero.

This function also contains two datasets:
1. Hi-C chromatin interaction matrix of chromosome 18. The data is based on mouse embryonic stem cell genome, and binning resolution of [40 kb](http://sysbio.rnet.missouri.edu/3dgenome/GSDB_old/details.php?id=TE1402WS) (Dixon *et al.* 2012).
2. Matched motif output by FIMO after scanning above chromsome 18 with 369 mouse motifs recently reported to regulate histone modification [Dataset_S01](https://www.pnas.org/content/116/9/3668/tab-figures-data) (Ngo *et al.*, 2019).

Finally, function `LinearizeNestedList` was written by [Akhil S. Bhel](https://github.com/mrdwab/mrdwabmisc/) (Mahto, Turner, and Bhel, 2016).

## References
Dixon JR, Selvaraj S, Yue F, Kim A *et al.* (2012) Topological domains in mammalian genomes identified by analysis of chromatin interactions. *Nature* **485**(7398):376--80. PMID: [22495300](https://www.ncbi.nlm.nih.gov/pubmed/22495300)

Ngo V, Chen Z, Zhang K, Whitaker J, Wang M, Wang W. (2019) Epigenomic analysis reveals DNA motifs regulating histone modifications in human and mouse. *PNAS* **116**(9): 3668--3677.

Mahto A, Turner R, Bhel AS. (2016) mrdwabmisc: Miscellaneous R functions, mostly for data processing. https://github.com/mrdwab/mrdwabmisc/.
