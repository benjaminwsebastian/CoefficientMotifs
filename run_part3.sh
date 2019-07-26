#!/bin/bash
# This script extracts the motifs from the meme output and translates them back to numbers
for i in {1..6}
do
    echo "Main coefficient motif for n = $i"
    cd "$i/memeOutput"
    python ../../scripts/translateMotif.py --f meme.html --bins ../"n$i"_"n100"_"bins.json"
    cd ../..
done
