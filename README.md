# Integrative genotype–environment association and single-cell transcriptomic analyses identify candidate mechanisms of environmental adaptation in pigs

This repository contains the analysis scripts for Li et al. (2026), published in *Molecular Biology and Evolution*, 43(10), msag260. [Journal article](https://academic.oup.com/mbe/article/43/10/msag260/8893486) · [DOI: 10.1093/molbev/msag260](https://doi.org/10.1093/molbev/msag260)

## Study overview

The study used whole-genome sequencing data from 1,151 pigs across 146 populations to investigate environmental adaptation. Three genotype–environment association (GEA) methods, redundancy analysis (RDA), BayPass, and latent factor mixed models (LFMM), were used to identify genomic regions associated with environmental variation.

Lung single-cell RNA sequencing of Tibetan and commercial crossbred pigs highlighted *ASXL1* as a candidate gene for high-altitude adaptation and implicated NOTCH signaling as a candidate pathway. The study also analyzed bulk lung RNA-seq data to compare candidate gene expression between high- and low-altitude groups.

## Repository structure

| Directory | Analyses |
| --- | --- |
| [Population analysis](./1.Population%20analysis/) | Neighbor-joining tree, principal component analysis, ancestry inference, linkage disequilibrium decay, runs of homozygosity, and selection scans |
| [Genotype–environment association](./2.GEAS/) | RDA, BayPass, and LFMM |
| [Single-cell RNA-seq analysis](./3.Single-cell%20RNA%20sequencing%20analysis/) | Preprocessing and integration, cell type annotation, trajectory inference, cell–cell communication, pseudobulk differential expression, enrichment, and in silico knockout analysis |
| [Bulk RNA-seq analysis](./4.RNAseq_analysis/) | Read alignment, duplicate marking, gene-level quantification, and expression visualization |

See the Methods section of the [journal article](https://academic.oup.com/mbe/article/43/10/msag260/8893486) for software versions, analysis parameters, and filtering criteria.

## Data availability

- **Lung single-cell RNA-seq data generated in this study** are deposited at CNCB-NGDC under BioProject [PRJCA058802](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA058802), which links to Genome Sequence Archive accession `CRA039139`.
- **Whole-genome sequencing and bulk lung RNA-seq sample information** is listed in Tables S1 and S12, respectively, in the supplementary material on the [journal article page](https://academic.oup.com/mbe/article/43/10/msag260/8893486).

## Citation

If you use this code, please cite the article.

> Li, D., Yuan, T., Li, X., Quan, Y., Tian, X., Cheng, Y., Wang, Y., Wu, X., Sun, B., Liu, W., Jiang, Y., Ding, R., Groenen, M. A. M., Derks, M. F. L., & Yu, T. (2026). Integrative genotype–environment association and single-cell transcriptomic analyses identify candidate mechanisms of environmental adaptation in pigs. *Molecular Biology and Evolution*, 43(10), msag260. [10.1093/molbev/msag260](https://doi.org/10.1093/molbev/msag260)

```bibtex
@article{li2026pigadaptation,
  author  = {Li, Dong and Yuan, Tiantian and Li, Xiang and Quan, Yu and
             Tian, Xuekai and Cheng, Ye and Wang, Yulong and Wu, Xinglong and
             Sun, Boyuan and Liu, Wansheng and Jiang, Yu and Ding, Rongrong and
             Groenen, Martien A. M. and Derks, Martijn F. L. and Yu, Taiyong},
  title   = {Integrative genotype--environment association and single-cell
             transcriptomic analyses identify candidate mechanisms of
             environmental adaptation in pigs},
  journal = {Molecular Biology and Evolution},
  year    = {2026},
  volume  = {43},
  number  = {10},
  eid     = {msag260},
  doi     = {10.1093/molbev/msag260}
}
```
