#!/usr/bin/bash

python fasta_to_mel_pairwise.py $1

for specie in $(cat species_list.txt); do
  ./muscle3.8.31_i86linux64 -in $specie -out ${specie::-6}_aligned.fasta -quiet 
done

python combine_alignments.py *pairwise_aligned.fasta

rm species_list.txt
rm *pairwise_aligned.fasta
rm *pairwise.fasta
