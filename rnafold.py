import pandas as pd
import RNA
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import argparse
import logging
from logging.handlers import QueueHandler, QueueListener
import gzip
import shutil
import os
import sys


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


def fold_fasta(batch, proc_id) -> pd.DataFrame:
    """ This function takes a batch of sequences and folds them using RNAFold, returning a dataframe with the sequence id, sequence, and predicted secondary structure"""
    solved_seqs = []
    for record in batch:
       (ss, mfe) = RNA.fold(str(record.seq))
       solved_seqs.append({'id':record.id, 'seq':record.seq, 'ss':ss, 'mfe':mfe})
       logging.info(f'Worker {proc_id} has folded a sequence')
    
    df = pd.DataFrame(solved_seqs)
    logging.info(f'Worker: {proc_id} has finished a {len(batch)} sequence batch and is exiting')
    return df
    
def worker_init(log_queue):
    qh = QueueHandler(log_queue)
    root = logging.getLogger()
    root.handlers = []
    root.addHandler(qh)
    root.setLevel(logging.INFO)

def main(args):
    n_workers = args.c
    file_path = args.f
    batches = []
    futures = []
    dfs = []
    log_queue = multiprocessing.Queue()
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s %(process)d %(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    listener = QueueListener(log_queue, handler)
    listener.start()

    if not os.path.exists(file_path):
        logging.critical('No Files were found at the specified path, this error is not recoverable, exiting now')
        raise FileNotFoundError

    if os.cpu_count() < args.c:
        logging.warning(f'CPU has less cores than requested, surplus processes will be queued, REQUESTED={args.c}, AVAIL={os.cpu_count()} \n')
    elif os.cpu_count() is None:
        logging.warning(f'Unable to determine number of CPU cores, proceeding with requested number of processes, REQUESTED={args.c} \n')
        
    with open(file_path) as f:
        
        for n, batch in enumerate(batch_iterator(SeqIO.parse(f, 'fasta'))):
            batches.append(batch)
        
    with ProcessPoolExecutor(max_workers=n_workers, initializer=worker_init, initargs=(log_queue, )) as scheduler:
        
        logging.info('Scheduling Jobs')
        n = 0
        for i in range(len(batches)):
            futures.append(scheduler.submit(fold_fasta, batches[i], n))
            logging.info(f'Queued Worker {n} and assigned {len(batches[i])} sequences to fold')
            n += 1

        for future in futures:
            dfs.append(future.result())

    listener.stop()

    dataset = pd.concat(dfs)

    if '.csv' in args.o:
        dataset.to_csv(args.o)
    elif '.tsv' in args.o:
        dataset.to_csv(args.o, sep='\t', index=False)
    elif '.json' in args.o:
        dataset.to_json(args.o)
    elif '.parquet' in args.o:
        dataset.to_parquet(args.o)


    if args.g and '.parquet' in args.o:
        logging.warning('Parquet output selected but gzip compression flag raised, gzip compression will be skipped')
        sys.exit(2)


    if args.g:
        
        with open(args.o, 'rb') as input_stream:
            with gzip.open(f'{args.o}.gz', 'wb') as output_stream:
                shutil.copyfileobj(input_stream, output_stream)
        
        if os.path.exists(args.o):
            os.remove(args.o)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='rnafold',
        description='a cli wrapper for Vienna RNAFold for RNA secondary structure predictions using MFE estimations'
    )

    parser.add_argument('-f', help='the directory containing fasta files', required=False, default='/home/FCAM/asands/projects/nonb/fasta/SGNex_Hct116_directRNA_replicate1_run1.fasta')
    parser.add_argument('-c', help='number of parallel processes to run folds with, default is 4', default=4, required=False, type=int)
    parser.add_argument('-o', help='the output file in txt format, supports .csv, .tsv, .json, and .parquet', default='output.csv')
    parser.add_argument('-v', help='determines the level of output from the script, adding this flag will print all thread logs', action='store_true')
    parser.add_argument('-g', help='optional flag to trigger gzip compression, this will not work if parquet format is already selected', action='store_true')

    arglist = parser.parse_args()

    if arglist.v:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    main(arglist)
