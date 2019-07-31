#!/bin/bash
# Create the data, then fit the polynomials, create figures, and fasta file
for i in {1..6}
do
    mkdir $i
    cd $i
    python3 ../scripts/create_data.py --n $(($i))
    python3 ../scripts/ls_fit.py --f "n$i.txt" --o "n$i"
    cd ..
done

