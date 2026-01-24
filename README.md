# nichenetr

<!-- badges: start -->
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/nichenetr)](https://zaoqu-liu.r-universe.dev/nichenetr)
[![R-CMD-check](https://github.com/Zaoqu-Liu/nichenetr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zaoqu-Liu/nichenetr/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://codecov.io/gh/Zaoqu-Liu/nichenetr/branch/master/graph/badge.svg)](https://codecov.io/gh/Zaoqu-Liu/nichenetr)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7074291.svg)](https://doi.org/10.5281/zenodo.7074291)
[![Nature Methods](https://img.shields.io/badge/Nature%20Methods-10.1038%2Fs41592--019--0667--5-blue)](https://www.nature.com/articles/s41592-019-0667-5)
<!-- badges: end -->

📖 **Documentation**: https://zaoqu-liu.github.io/nichenetr/

## Overview

**NicheNet** is a computational framework for modeling intercellular communication by linking extracellular ligands to downstream intracellular gene regulatory changes. Unlike conventional ligand-receptor inference methods that solely identify potential receptor-ligand interactions based on co-expression patterns, NicheNet integrates prior knowledge of signaling pathways and gene regulatory networks to predict which ligands from sender cells are most likely to affect gene expression in receiver cells.

### Key Publications

> Browaeys R, Saelens W, Saeys Y. **NicheNet: modeling intercellular communication by linking ligands to target genes.** *Nature Methods* 17, 159-162 (2020). [doi:10.1038/s41592-019-0667-5](https://doi.org/10.1038/s41592-019-0667-5)

> **Unraveling cell-cell communication with NicheNet by inferring active ligands from transcriptomics data.** *Nature Protocols* (2024). [doi:10.1038/s41596-024-01121-9](https://doi.org/10.1038/s41596-024-01121-9)

## Installation

### From R-universe (Recommended)

```r
install.packages("nichenetr", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("Zaoqu-Liu/nichenetr")
```

### System Requirements

- **R version**: >= 4.0.0
- **Operating systems**: Windows, macOS, Linux
- **Dependencies**: Seurat (v4/v5), tidyverse, igraph, Matrix, Rcpp

## Quick Start

```r
library(nichenetr)
library(Seurat)
library(tidyverse)

# Load pre-built prior model
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

# Run NicheNet analysis
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seurat_obj,
  receiver = "Receiver_CellType",
  sender = c("Sender_CellType1", "Sender_CellType2"),
  condition_oi = "Treatment",
  condition_reference = "Control",
  condition_colname = "condition",
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)
```

## Tutorials

| Vignette | Description |
|----------|-------------|
| [Algorithm Overview](https://zaoqu-liu.github.io/nichenetr/articles/algorithm_overview.html) | Comprehensive overview of the NicheNet algorithm |
| [Step-by-step Analysis](https://zaoqu-liu.github.io/nichenetr/articles/seurat_steps.html) | Detailed walkthrough with Seurat |
| [Seurat Wrapper](https://zaoqu-liu.github.io/nichenetr/articles/seurat_wrapper.html) | Streamlined analysis using wrapper functions |

See the [full documentation](https://zaoqu-liu.github.io/nichenetr/) for more vignettes.

## Data Resources

Pre-built models available at [Zenodo](https://zenodo.org/record/7074291):

- Ligand-Target Matrix (Human/Mouse)
- Weighted Networks (Human/Mouse)
- Ligand-Receptor Networks (Human/Mouse)

## Citation

```bibtex
@article{browaeys2020nichenet,
  title={NicheNet: modeling intercellular communication by linking ligands to target genes},
  author={Browaeys, Robin and Saelens, Wouter and Saeys, Yvan},
  journal={Nature Methods},
  volume={17},
  pages={159--162},
  year={2020},
  doi={10.1038/s41592-019-0667-5}
}
```

## License

GPL-3.0

## Contact

For questions, please open an issue on [GitHub](https://github.com/Zaoqu-Liu/nichenetr/issues).
