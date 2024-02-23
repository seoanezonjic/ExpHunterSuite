
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExpHunterSuite

<!-- badges: start -->

[![R-CMD-check](https://github.com/seoanezonjic/ExpHunterSuite/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seoanezonjic/ExpHunterSuite/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## DOCUMENTATION - VIGNETTES

[Link to ExpHunterSuite Vignettes](https://bioconductor.org/packages/3.18/workflows/vignettes/ExpHunterSuite/inst/doc/hunter.html)

## DOCUMENTATION - README

# DEgenes Hunter

> Protocol to undertake the key steps in the analysis of transcriptomics
> experiments, including differentially expressed gene detection, a
> co-expression analysis and the functional enrichment.

## REQUIREMENTS

  - R version 4

## INSTALLATION

### Install latest R-Version

Go to page <https://cran.r-project.org> and install the latest R version
on your computer.

### Install DEgenes Hunter from console (using devtools)

To install DEgenes Hunter from console, we use the devtools utility to
install R packages from GitHub. Please type the following
commands:

``` bash
    Rscript -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'
    Rscript -e 'devtools::install_github("seoanezonjic/ExpHunterSuite", dependencies=TRUE)'
```

Sometimes reactome.db download fails because this package is too large
(\>400 Mb) and R has a download timeout of 60 seconds by default.

``` bash
    Rscript -e "getOption('timeout')"
        [1] 60
```

This error can be solved by installing reactome.db from source or
extending the timeout threshold before installing the ExpHunterSuite
package with the
command:

``` bash
    Rscript -e 'options(timeout=1000); devtools::install_github("seoanezonjic/ExpHunterSuite")'
```

Then, please create the folder where you want to install the command
line scripts, copy the scripts there and make them command line
accesible using these commands:

``` bash
    mkdir install_folder
    Rscript -e "ExpHunterSuite::install_DEgenes_hunter('install_folder')"
    export PATH=path_to_install_folder:$PATH
```

This export PATH can also be added to the .bashrc or .bash\_profile
files.

## REFERENCES AND CITATION

Please cite ExpHunterSuite as:



Jabato, F., Córdoba-Caballero, J., Rojano, E., Romá-Mateo, C., Sanz, P., Pérez, B., Gallego, D., Seoane, P., Ranea, J., and Perkins, J. 2021. Gene expression analysis method integration and co-expression module detection applied to rare glucide metabolism disorders using ExpHunterSuite. *Scientific Reports* 2021 11:1, 11(1), p.1–12. <doi:http://dx.doi.org/0.1038/s41598-021-94343-w>


Isabel González Gayte, Rocío Bautista Moreno, Pedro Seoane Zonjic and Manuel G. Claros (2017). DEgenes Hunter - A Flexible R Pipeline for Automated RNA-seq Studies in Organisms without Reference Genome. *Genomics And Computational Biology* , 3(3), e31. <doi:http://dx.doi.org/10.18547/gcb.2017.vol3.iss3.e31>


