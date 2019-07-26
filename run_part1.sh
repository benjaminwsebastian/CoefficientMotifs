#!/bin/bash
# Create the data, then fit the polynomials, create figures, and fasta file
for i in {1..6}
do
    mkdir $i
    cd $i
    python ../scripts/createData.py --n $(($i))
    python ../scripts/lsFit.py --f "n$i.txt" --o "n$i"
    cd ..
done

