#!/bin/bash
for i in {1..6}
do
    mkdir $i
    cd $i
    python -W ignore ../scripts/createData.py --n $(($i))
    python -W ignore ../scripts/lsFit.py --f "n$i.txt" --o "n$i"
    cd ..
done

