#!/bin/bash

cd Data
ssu-align barrnap_longest_16S.fasta myseqs
ssu-mask myseqs
ssu-mask --stk2afa myseqs
ssu-draw myseqs
#convert masked stk file to fasta
cd myseqs
ssu-mask -a --stk2afa myseqs.bacteria.mask.stk