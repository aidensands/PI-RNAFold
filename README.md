# PI-RNAFold
This is a parallel implementation of ViennaRNA's RNAFold. Since the existing tool is CPU bound, this script dramatically speeds up compute time using multiple CPU cores. All credit goes to the creators of ViennaRNA.

## New
I'm currently migrating this codebase to Rust! The python section of the code is still perfectly functional and recently recieved memory performance improvements. If you wish to use the python version still or the rust version is not
ready yet, all you have to do is cd into the python folder or run the script using:
```
python legacy/rnafold.py --help
```


## Installation
To remotely download the script ensure that you have Python installed and the use:
```
git clone https://github.com/aidensands/PI-RNAFold.git
```
Then, install the scripts requirements using 
```
python -m venv [env name]
pip install -r requirements.txt
```
## Usage
To use this tool as intended you must have your sequences in .fasta/.fa format. The command to run PI-RNAFold is as follows:
```
python rnafold.py \
    -f [Path to your fasta file]
    -o [Name for output file in supported format (.csv, .tsv, .json, .parquet)]
    -c [Number of consecutive processes to handle folding]
    -v [Raise this flag for verbose output (no parameter needed)]
    -g [Raise this flag to trigger gzip compression of output (no parameter needed)]
```

## Effeciency
This script avoids the dangers of concurrent processing by using Python's pre-existing concurrency libraries. More specifically all processes are handled by a ProcessPoolExecuter scheduler that prevents us from having to worry about shared memory race conditions. The official implementation of ViennaRNA's RNAFold is:

$$O(N \cdot L^3)$$
$$N:\text{ Number of Sequences}$$
$$L: \text{ Length of Sequence}$$

The new parallel implementation's speedup is dependant on your resources, making it especially useful for teams working in an HPC environment. The new complexity is embarassingly parallel and thus its complexity is listed below:

$$O(\frac{N\cdot L^3}{P})$$
$$N:\text{ Number of Sequences}$$
$$L: \text{ Length of Sequence}$$
$$P: \text{ Number of Processor Cores}$$