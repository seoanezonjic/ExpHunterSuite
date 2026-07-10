#! /usr/bin/env bash

if [[ "$DEGHUNTER_MODE" == "" ]]; then
	. ~soft_bio_267/initializes/init_degenes_hunter
fi

download_KEGG_file.R -O Human
download_KEGG_file.R -O Mouse
