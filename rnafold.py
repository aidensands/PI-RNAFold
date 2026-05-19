import pandas as pd
import RNA
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import glob
import argparse
import logging
import os


def batch_iterator(iterator, batch_size = 100):
    """ This function builds a generator that can be iterated through
        it contains sequence batches of size 'batch_size'
    """
    batch = list()

    for entry in iterator:
    
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    
    if batch:
        yield batch


def fold_fasta(batch, proc_id):
    """ This function takes a batch of sequences and folds them using RNAFold, returning a dataframe with the sequence id, sequence, and predicted secondary structure"""
    solved_seqs = []
    for record in batch:
       (ss, mfe) = RNA.fold(str(record.seq))
       solved_seqs.append({'id':record.id, 'seq':record.seq, 'ss':ss, 'mfe':mfe})
       logging.info(f'Worker {proc_id} has folded a sequence')
    
    df = pd.DataFrame(solved_seqs)
    logging.info(f'Worker: {proc_id} has finished a {len(batch)} sequence batch and is exiting')
    return df
    

def main(args):
    n_workers = args.c
    file_glob = glob.glob(args.d + '*.fasta')
    batches = []
    futures = []
    dfs = []

    if not file_glob:
        logging.critical('No Files were found at the specified path, this error is not recoverable, exiting now')
        raise FileNotFoundError

    if os.cpu_count() < args.c:
        logging.warning(f'CPU has less cores than requested, surplus processes will be queued, REQUESTED={args.c}, AVAIL={os.cpu_count()} \n')

    for file in file_glob:
        
        with open(file) as f:
            
            for n, batch in enumerate(batch_iterator(SeqIO.parse(open(file), 'fasta'))):
                batches.append(batch)
            
            with ProcessPoolExecutor(max_workers=n_workers) as scheduler:
                logging.info('Scheduling Jobs')
                n = 0
                for i in range(len(batches)):
                    futures.append(scheduler.submit(fold_fasta, batches[i], n))
                    logging.info(f'Queued Worker {n} and assigned {len(batches[i])} sequences to fold')
                    n += 1
                    
                for future in futures:
                    dfs.append(future.result())
        
            f.close()
        
        batches.clear()

        
    dataset = pd.concat(dfs)
    dataset.to_csv(args.o)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='rnafold',
        description='a cli wrapper for Vienna RNAFold for RNA secondary structure predictions using MFE estimations'
    )

    parser.add_argument('-d', help='the directory containing fasta files', required=False, default='/home/FCAM/asands/projects/nonb/fasta/SGNex_Hct116_directRNA_replicate1_run1.fasta')
    parser.add_argument('-c', help='number of parallel processes to run folds with, default is 4', default=4, required=False, type=int)
    parser.add_argument('-o', help='the output file in txt format', default='output.csv')
    parser.add_argument('-v', help='determines the level of output from the script, adding this flag will print all thread logs', action='store_true')

    arglist = parser.parse_args()

    if arglist.v:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    main(arglist)
