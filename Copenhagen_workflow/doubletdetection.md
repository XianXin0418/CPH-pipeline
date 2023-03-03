Doublet detection using DoubletDetection
================
Xian Xin
2023-03-02

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#load-count-matrix" id="toc-load-count-matrix">Load count
    matrix</a>
-   <a href="#default-pipeline" id="toc-default-pipeline">Default
    pipeline</a>
-   <a href="#export-doubletdetection-results-to-r"
    id="toc-export-doubletdetection-results-to-r">Export DoubletDetection
    results to R</a>

## Introduction

[DoubletDetection](https://github.com/JonathanShor/DoubletDetection) can
identify doublets in scRNA-seq data. Similar as *scrublet*, It requires
count matrix (cells as rows and genes as columns) as input to calculate
doublet scores and call doublets in the sample. It is recommended to
**run on each sample separately** but not on merged count matrix, and
the **dataset should includes several cell types** (optional but highly
recommended). After finishing doublet detection, we also recommend to
check and adjust (if necessary) the doublet score threshold.

In this vignette, we use the *CellBender* filtered count matrix from one
SCN2A mouse sample to demonstrate how to run *DoubletDetection* in
reticulate Python environment. This vignette is based on
[*DoubletDetection*
tutorial](https://doubletdetection.readthedocs.io/en/latest/tutorial.html).

For detailed instruction of *DoubletDetection* installation, please
refer to [Quality Control using CRMetrics vignette](./QC.md).

## Load count matrix

To run Python in R Markdown, we need to use *reticulate* to load the
conda environment which contains *DoubletDetection* module.

``` r
library(magrittr)
library(reticulate)
use_miniconda(condaenv = "/home/gjl413/.conda/envs/r-crmetrics/bin/python")
```

To load the count matrix, we simply use a *Seurat* function to read the
hdf5 file.

``` r
cm <- Seurat::Read10X_h5(filename = "/people/gjl413/data/FORpipeline_example/P14/P14_het_1/outs/cellbender_filtered.h5", use.names = TRUE)
```

## Default pipeline

**In the rest part of this vignette, the codes are all in Python**. For
running Python in R Markdown, please refer to [*reticulate*
document](https://rstudio.github.io/reticulate/).

``` python
import numpy as np
import doubletdetection
import matplotlib.pyplot as plt
```

Run *DoubletDetection* with default parameters except `n_jobs`.

We can easily use objects from R environment in Python with `r.VAR_NAME`
(and use Python objects in R with `py$VAR_NAME`).

``` python
clf = doubletdetection.BoostClassifier(
  clustering_algorithm='phenograph', # One of [“louvain”, “leiden” "phenograph"]. "louvain" and leiden refer to the scanpy module.
  n_iters=10, # Number of fit operations from which to collect p-values
  pseudocount=0.1, # Pseudocount used in normalize_counts. If 1 is used, and standard_scaling=False, the classifier is much more memory efficient; however, this may result in fewer doublets detected.
  standard_scaling=False, # Set to True to enable standard scaling of normalized count matrix prior to PCA. Recommended when not using Phenograph. Defaults to False.
  n_jobs=32 # Number of jobs to use. Speeds up neighbor computation.
)
```

`doublet_labels` is a one-dimensional `numpy ndarray` with the value `1`
representing a detected doublet, `0` a singlet, and `np.nan` an
ambiguous cell.

``` python
doublet_labels = clf.fit(r.cm).predict(
  p_thresh=1e-07, # Hypergeometric test p-value threshold that determines per iteration doublet calls.
  voter_thresh=0.9 # Fraction of iterations a cell must be called a doublet
)
```

    ##   0%|          | 0/10 [00:00<?, ?it/s] 10%|#         | 1/10 [01:37<14:41, 97.91s/it] 20%|##        | 2/10 [03:50<15:46, 118.29s/it] 30%|###       | 3/10 [05:56<14:13, 121.95s/it] 40%|####      | 4/10 [08:29<13:24, 134.17s/it] 50%|#####     | 5/10 [10:27<10:41, 128.22s/it] 60%|######    | 6/10 [12:18<08:09, 122.36s/it] 70%|#######   | 7/10 [14:05<05:51, 117.29s/it] 80%|########  | 8/10 [15:58<03:52, 116.08s/it] 90%|######### | 9/10 [17:49<01:54, 114.36s/it]100%|##########| 10/10 [19:52<00:00, 116.96s/it]100%|##########| 10/10 [19:52<00:00, 119.20s/it]

`doublet_scores` is a one-dimensional numpy ndarray representing a score
for how likely a cell is to be a doublet. The score is used to create
the labels.

``` python
doublet_scores = clf.doublet_score()
```

## Export DoubletDetection results to R

``` r
# execute in R
doublet.scores <- py$doublet_scores %>% setNames(colnames(cm))
doublet.labels <- py$doublet_labels %>% setNames(colnames(cm))
```

``` r
head(doublet.scores)
```

    ## TTGCCTGCAGCGTGAA-1 TGCTCCAAGCCACTCG-1 CATCCCAGTCCCTAAA-1 TGCATGAAGGTTGTTC-1 
    ##          50.027771           0.000000          54.053625           0.000000 
    ## TGGAACTAGGTCCGAA-1 TGCGACGCATAACTCG-1 
    ##           1.960648           2.049571

``` r
head(doublet.labels)
```

    ## TTGCCTGCAGCGTGAA-1 TGCTCCAAGCCACTCG-1 CATCCCAGTCCCTAAA-1 TGCATGAAGGTTGTTC-1 
    ##                  1                NaN                  1                NaN 
    ## TGGAACTAGGTCCGAA-1 TGCGACGCATAACTCG-1 
    ##                  0                  0

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux 8.7 (Ootpa)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/local/R-4.2.2/lib64/R/lib/libRblas.so
    ## LAPACK: /usr/local/R-4.2.2/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] reticulate_1.28 magrittr_2.0.3 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Seurat_4.3.0           Rtsne_0.16             colorspace_2.1-0      
    ##   [4] deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.4        
    ##   [7] rprojroot_2.0.3        rstudioapi_0.14        spatstat.data_3.0-0   
    ##  [10] leiden_0.4.3           listenv_0.9.0          bit64_4.0.5           
    ##  [13] ggrepel_0.9.2          fansi_1.0.4            codetools_0.2-18      
    ##  [16] splines_4.2.2          knitr_1.42             polyclip_1.10-4       
    ##  [19] jsonlite_1.8.4         ica_1.0-3              cluster_2.1.4         
    ##  [22] png_0.1-8              uwot_0.1.14            shiny_1.7.4           
    ##  [25] sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.2.2        
    ##  [28] httr_1.4.4             SeuratObject_4.1.3     Matrix_1.5-3          
    ##  [31] fastmap_1.1.0          lazyeval_0.2.2         cli_3.6.0             
    ##  [34] later_1.3.0            htmltools_0.5.4        tools_4.2.2           
    ##  [37] igraph_1.3.5           gtable_0.3.1           glue_1.6.2            
    ##  [40] RANN_2.6.1             reshape2_1.4.4         dplyr_1.1.0           
    ##  [43] rappdirs_0.3.3         Rcpp_1.0.10            scattermore_0.8       
    ##  [46] vctrs_0.5.2            nlme_3.1-160           spatstat.explore_3.0-6
    ##  [49] progressr_0.13.0       lmtest_0.9-40          spatstat.random_3.1-3 
    ##  [52] xfun_0.37              stringr_1.5.0          globals_0.16.2        
    ##  [55] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
    ##  [58] irlba_2.3.5.1          goftest_1.2-3          future_1.31.0         
    ##  [61] MASS_7.3-58.1          zoo_1.8-11             scales_1.2.1          
    ##  [64] promises_1.2.0.1       spatstat.utils_3.0-1   parallel_4.2.2        
    ##  [67] RColorBrewer_1.1-3     yaml_2.3.7             pbapply_1.7-0         
    ##  [70] gridExtra_2.3          ggplot2_3.4.0          stringi_1.7.12        
    ##  [73] rlang_1.0.6            pkgconfig_2.0.3        matrixStats_0.63.0    
    ##  [76] evaluate_0.20          lattice_0.20-45        ROCR_1.0-11           
    ##  [79] purrr_1.0.1            tensor_1.5             patchwork_1.1.2       
    ##  [82] htmlwidgets_1.6.1      bit_4.0.5              cowplot_1.1.1         
    ##  [85] tidyselect_1.2.0       here_1.0.1             parallelly_1.34.0     
    ##  [88] RcppAnnoy_0.0.20       plyr_1.8.8             R6_2.5.1              
    ##  [91] generics_0.1.3         DBI_1.1.3              pillar_1.8.1          
    ##  [94] withr_2.5.0            fitdistrplus_1.1-8     survival_3.4-0        
    ##  [97] abind_1.4-5            sp_1.6-0               tibble_3.1.8          
    ## [100] future.apply_1.10.0    hdf5r_1.3.8            KernSmooth_2.23-20    
    ## [103] utf8_1.2.3             spatstat.geom_3.0-6    plotly_4.10.1         
    ## [106] rmarkdown_2.20         grid_4.2.2             data.table_1.14.6     
    ## [109] digest_0.6.31          xtable_1.8-4           tidyr_1.3.0           
    ## [112] httpuv_1.6.8           munsell_0.5.0          viridisLite_0.4.1
