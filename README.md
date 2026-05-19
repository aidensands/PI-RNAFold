# PI-RNAFold
This is a parallel implementation of ViennaRNA's RNAFold. Since the existing tool is CPU bound, this script dramatically speeds up compute time using multiple CPU cores. All credit goes to the creators of ViennaRNA 

## Usage
To use this tool as intended you must have your sequences in .fasta/.fa format. The command to run PI-RNAFold is as follows:
```
python rnafold.py \
    -f [Path to your fasta file]
    -o [Name for output file in csv format]
    -c [Number of consecutive processes to handle folding]
    -v [Raise this flag for verbose output (no parameter needed)]
```

