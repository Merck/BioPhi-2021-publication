#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import os
import numpy as np
from collections import OrderedDict
import gzip
from tqdm import tqdm
from abnumber import Chain
import re
import requests
import time
import sys

T20_REGEX = re.compile('<td>T20 Score:</td><td>([0-9.]+)</td>')

def get_t20_online(seq):
    chain = Chain(seq, scheme='imgt')
    chain_type = 'vh' if chain.chain_type == 'H' else ('vl' if chain.chain_type == 'L' else 'vk')
    html = None
    for retry in range(5):
        url = f'https://dm.lakepharma.com/cgi-bin/blast.py?chain={chain_type}&region=1&output=3&seqs={seq}'
        try:
            request = requests.get(url)
            if request.ok:
                html = request.text
                break
        except Exception as e:
            print(e)
        time.sleep(0.5 + retry * 5)
        print('Retry', retry+1)
    if not html:
        sys.exit(1)
    matches = T20_REGEX.findall(html)
    time.sleep(1)
    if not matches:
        print(html)
        raise ValueError(f'Error calling url {url}')
    return float(matches[0])

def get_multiple_t20_online(queries):
    results = []
    for query in tqdm(queries):
        t20 = get_t20_online(query.seq)
        results.append(OrderedDict(
            id=query.id,
            description=query.description,
            t20=t20
        ))
    return pd.DataFrame(results)

def get_multiple_t20_custom_db(queries, targets, num=20):
    from skbio.alignment import StripedSmithWaterman
    from skbio.alignment._pairwise import blosum50
    results = []
    for query in tqdm(queries):
        aligner = StripedSmithWaterman(
            str(query.seq), 
            protein=True,
            substitution_matrix=blosum50
        )
        best_identities = []
        best_aligned_query = None
        best_aligned_target = None
        
        for target in targets:
            alignment = aligner(target)
            query_length = len(alignment.query_sequence)
            identity = sum([a == b for a, b in zip(alignment.aligned_query_sequence, alignment.aligned_target_sequence)]) / query_length
            if len(best_identities) < num:
                best_identities.append(identity)
                best_identities = sorted(best_identities)
            elif identity > best_identities[0]:
                best_identities = best_identities[1:] + [identity]
                best_identities = sorted(best_identities)
            
            if identity == best_identities[-1]:
                best_aligned_query = alignment.aligned_query_sequence
                best_aligned_target = alignment.aligned_target_sequence
        result = OrderedDict(
            query_seq=best_aligned_query,
            id=query.id,
            description=query.description,
            best_seq=best_aligned_target,
            best_identity=best_identities[-1]
        )
        result[f'T{num}'] = np.mean(best_identities)
        results.append(result)
    return pd.DataFrame(results)
    
if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="Query FASTA sequence file.")
    parser.add_argument("output", help="Output TSV file path.")
    parser.add_argument("--db", help="Custom human gzipped TXT sequence database file.")
    parser.add_argument("-n", "--num", default=20, type=int, help="Number of closest sequences to evaluate.")

    options = parser.parse_args()
    
    
    queries = SeqIO.parse(options.query, 'fasta')
    
    if options.db:
        with gzip.open(options.target, 'rt') as f:
            targets = [l.strip() for l in f.readlines()]

        df = get_multiple_t20_custom_db(queries, targets, num=options.num)
    else:
        assert options.num == 20, 'You need to provide a custom --db to run with n != 20'
        print('Note: The sequences will be processed through lakepharma T20 service! Sleeping for 10s, press Ctrl+C to cancel...')
        time.sleep(10)
        print('Processing...')
        df = get_multiple_t20_online(queries)
    
    df.to_csv(options.output, sep='\t', index=False)
    print(f'Saved to: {options.output}')
