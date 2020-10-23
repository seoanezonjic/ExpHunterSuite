#! /usr/bin/env bash

. ~soft_bio_267/initializes/init_degenes_hunter

export DEGHUNTER_MODE=DEVELOPMENT
export PATH=../scripts:$PATH

which degenes_Hunter.R
# grep -o '"\-\-[A-Za-z0-9_]*' options_file | sed 's/"--\([A-Za-z0-9_]*\)/\1=opt$\1,/g'
# grep -o '"\-\-[A-Za-z0-9_]*' degenes_Hunter.R | sed 's/"--\([A-Za-z0-9_]*\)/#'\'' @param \1/g' To convert params to Roxygen doc
degenes_Hunter.R -h
