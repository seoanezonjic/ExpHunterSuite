#! /usr/bin/env bash

. ~soft_bio_267/initializes/init_degenes_hunter
# https://laderast.github.io/2019/02/12/package-building-description-namespace/
# Rscript -e "usethis::use_package('optparse')" # To add dependencies that are not consigned by roxygen documentation

current=`pwd`
cd ../../
Rscript -e "devtools::document()"
Rscript -e "devtools::install(dependencies=FALSE)"

cd $current
mkdir install_folder
Rscript -e "DEgenesHunter::install_DEgenes_hunter('install_folder')"
install_folder/degenes_Hunter.R -h
