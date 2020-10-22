#! /usr/bin/env bash

. ~soft_bio_267/initializes/init_degenes_hunter

export DEGHUNTER_MODE=DEVELOPMENT
export PATH=../scripts:$PATH

which degenes_Hunter.R
# grep -o '"\-\-[A-Za-z0-9_]*' options_file | sed 's/"--\([A-Za-z0-9_]*\)/\1=opt$\1,/g'
degenes_Hunter.R -h
