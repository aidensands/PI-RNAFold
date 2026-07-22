# PI-RNAFold
This is a parallel implementation of ViennaRNA's RNAFold. Since the existing tool is CPU bound, this script dramatically speeds up compute time using multiple CPU cores. All credit goes to the creators of ViennaRNA.

## New
This codebase has just been migrated to Rust! PI-RNAFold has just recieved massive memory improvements and performance speedups from this migration. The python section of the code is still perfectly functional and recently recieved memory performance improvements. If you wish to use the python version, all you have to do is cd into the python folder or run the script using:
```
python legacy/rnafold.py --help
```

## Prerequisites:
In order to use the rust version of PI-RNAFold correctly you must have the following prerequisites installed:
- Rust
- Cargo (natively installed with Rust)
- The viennarna package, this can be done using any package manager or through an environment manager like conda by installing it into your active environment.

## Installation and Usage:
Clone this repository using:
```
git clone https://github.com/aidensands/PI-RNAFold.git
```

Compile the PI-RNAFold binary:

```
cargo build --release
```
The code will compile and an executable will be produced under `target/release/`. From here you are ready to go! Try running:

```
parallel_rnafold --help
```

## Legacy

```
python -m venv [env name]
pip install -r requirements.txt
```
## Legacy Usage
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
