Single-cell RNA-seq data analysis pipeline
================
Xian Xin
2023-03-02

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#content" id="toc-content">Content</a>
-   <a href="#will-be-added" id="toc-will-be-added">Will be added</a>
    -   <a href="#tools-and-analyses" id="toc-tools-and-analyses">Tools and
        analyses</a>
    -   <a href="#other-works" id="toc-other-works">Other works</a>

## Introduction

This vignettes collection introduces the tools what we mainly use for
scRNA-seq analysis in KKH lab, including introduction of each tool,
installation guide, usage, output format and corresponding graph output
exhibition.

## Content

[1. Quality Control using CRMetrics](./Copenhagen_workflow/QC.md)

-   [Ambient RNA removal using
    CellBender](./Copenhagen_workflow/cellbender.md)

-   [Ambient RNA removal using SoupX](./Copenhagen_workflow/soupx.md)

-   [Doublet detection using
    scrublet](./Copenhagen_workflow/scrublet.md)

-   [Doublet detection using
    DoubletDetection](./Copenhagen_workflow/doubletdetection.md)

[2. Single-cell RNA-seq datasets processing using Pagoda2 and
Conos](./Copenhagen_workflow/align_cluster.md)

-   [Pagoda2 analysis of single
    sample](./Copenhagen_workflow/pagoda2.md)

[3. Comparative analysis between groups of scRNA-seq samples using
Cacoa](./Copenhagen_workflow/cacoa.md)

[4. Transcriptional factor activity and pathway activity inference using
decoupleR](./Copenhagen_workflow/decoupler.md)

## Will be added

### Tools and analyses

**Python implementation of decoupleR**: Automated cell type annotation,
pathway activity inference, TF activity inference, functional enrichment
of biological terms, pseudo-bulk functional analysis and spatial
functional analysis. Faster and memory efficient comparing running
decoupleR in R.

**slingshot and tradeSeq**: Trajectory inference and trajectory-based
differential expression analysis.

**CellRanger**: 10X reads mapping and gene expression quantification.

**CELLECT**: Integrating scRNA-seq data with genetics data.

**Signac**: Analysis of scATAC-seq data and multimodal integration.

**CellChat or NeuronChat**: Cell-cell communication inference – need
further test for *NeuronChat*.

**IReNA or SCENICplus**: GRN analysis integrating scRNA-seq and
scATAC-seq data – need further test for *SCENICplus*

### Other works

Need to add the vignette guiding data format conversion between
`Seurat`, `Conos/Pagoda2`, `AnnData` (for scanpy),
`SingleCellExperiment`, etc. For conversion between `Seurat`, `AnnData`
and `SingleCellExperiment`, you can refer to the [Satija Lab
vignette](https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html).

Need to add comparison among popular tools for each analysis (pros and
cons).

Add reference for each vignette. Refine format and English writing.
