#!/bin/bash
#define par regions								
#NOTE: THIS IS IN HG19!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par1_start=60001; par1_end=2699520; par2_start=154931044; par2_end=155260560;

#add par or nonpar based on location																
less -S Files/Tukiainen_xinact.tsv |  sed 's/\r//g' | awk -F"\t" -v par1_s=$par1_start -v par1_e=$par1_end -v par2_s=$par2_start -v par2_e=$par2_end \
'{if($4>=par1_s && $5 <=par1_e) print $0"\tPAR1\tPAR"; else if($4>=par2_s && $5<=par2_e) print $0"\tPAR2\tPAR"; else if ($1=="Gene name"){print $0"\tPAR_STATUS\tPAR_BINARY"}else print $0"\tNONPAR\tNONPAR"}' \
> Files/Tukiainen_xinact_par.tsv

echo "DONE"


