# Spatial heterogeneity of tumor microenvironment influences the prognosis and therapeutic response of clear cell renal cell carcinoma patients

# Repo Contents
- [R](./R): `R` package code.
- [example_data](./example_data): test data
- [1Cellannotate](./1Cellannotate): result of cell cluster and annotation 
- [2Mixscore](./2Mixscore): mix score of the sample
- [3CellNeighborhood](./3CellNeighborhood): result of cell neighborhood 
## Install software and R package

### Installing R version 4.0.5 on CentOS 7.9.2009

the latest version of R can be installed by adding the latest repository to `yum`:

```
tar -xzvf R-4.0.5.tar.gz
yum group install "Development tools"
yum install readline-devel
yum install xz xz-devel 
yum install pcre pcre-devel
yum install libcurl-devel
yum install texlive
yum install java-1.8.0-openjdk
yum install *gfortran*
yum install zlib*
yum install bzip2-*
chown -R root:root R-4.0.5/
./configure –with-x=no
```

### Package Installation
Users should install the following packages 

```
install.packages(c('Seurat', 'RColorBrewer', 'ggthemes', 'ggplot2', 'viridis', 'S4Vectors', 'SingleCellExperiment', 'SpatialExperiment','magrittr','dplyr'))
```
The versions of software are, specifically:
```
seurat 4.0.5
ggplot2 3.3.5
RColorBrewer 1.1-2
ggthemes 4.2-4
viridis 0.6.2
S4Vectors 0.28.1
SingleCellExperiment 1.12.0
SpatialExperimentc 1.0.0
magrittr 2.0.1
dplyr 2.1.1
```

#### run example data

Firstly, cell cluster dimension reduction and cell identity was defined by seurat R package.
Secondly,caculate mixscore using mixing_score_summary function by SPIAT R package
Finally,spatial data was used to identify cell neighborhoods

```
Rscript scripts/IMC.step1.annotateCell.R example_data/P0_histocat_outfile.csv example_data/panel.csv 1Cellannotate/ P0
Rscript scripts/IMC.step2.mixscore.R example_data/P0_histocat_outfile.csv 1Cellannotate/P0.sc.cellAnnotate.Rds 2Mixscore/ P0
Rscript scripts/IMC.step3.identifyCN.R example_data/P0_histocat_outfile.csv 1Cellannotate/P0.sc.cellAnnotate.Rds 2Mixscore/P0.mixscore.csv example_data/panel.csv 3CellNeighborhood/ P0
```
### Input File Description
- *_histocat_outfile.csv：CSV output file of histocat software
- panel.csv: panel data 
- *.sc.cellAnnotate.Rds: output file after running IMC.step1.annotateCell.R
- *.mixscore.csv: output file after running IMC.step2.mixscore.R

Note: SPIAT R package(https://bioconductor.org/packages/release/bioc/html/SPIAT.html) source code was download under the script/SPIAT-main/R 
directory, this package was called on our script
