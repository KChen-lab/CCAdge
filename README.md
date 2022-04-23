**HighOrder_SComparison** is an R package for detecting differenital genes associated with phenotype from single cell profiles. This tool is developed and maintained by [Ken chen's lab](https://www.mdanderson.org/research/departments-labs-institutes/labs/ken-chen-laboratory.html) in MDACC. There is a pressing need for computatinal tools to enable the detection of genes associated with phenotypes in single-cell level. This task is very challening in single-cell studies due to the existence of: 1) changes in cell type composition among samples and 2) unknow technological batch effects. The proposed high-order gene detection method first leverages the canonical correlation analysis (CCA) to align cell populations from two sample into shared space, and then detects gene changed in each latent space. The gene difference in latent spaces are finally aggregated into the population-level mismatching score, which can capture gene change due to subtle cell states.

<image src="./doc/image/logo.png" width="800"> 
  
  
## Version History 
### v1.0.0 [4/22/2022]
* Release test version.
  
  


## System Requirements

### Hardware requirements
The package requires only a standard computer with enough RAM to support the in-memory operations. For minimal performance, please make sure that the computer has at least about `10 GB` of RAM. For optimal performance, we recommend a computer with the following specs:

* RAM: 10+ GB
* CPU: 4+ cores, 2.3 GHz/core
* Seurat (>=V3.0)
* irlba, umap, rdist, EnhancedVolcano, gtools, ggpubr

## Usage 

For usage examples and guided walkthroughs, check the `vignettes` directory of the repo.

*  [Quick start in spatial transcriptomics data](https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/ST_DEG.html)

*  More to be added ... 


## Bug report

## License
This project is covered under the **GNU General Public License 3.0**.

## Citation


