# Pipe a cat'd file into this! eg.
# cat myResidues.residue | sortResidues.sh
sort -t: -k 1,2 -k 3,3n
