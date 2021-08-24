#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import os
import numpy as np
from collections import OrderedDict
from tqdm import tqdm
from abnumber import Chain
import re
import requests
import time

SCORE_REGEX = re.compile('<h3>The Z-score value of the Query sequence is: (-?[0-9.]+)</h3>')

def get_z_score_online(seq):
    chain = Chain(seq, scheme='imgt')
    chain_type = 'human_heavy' if chain.chain_type == 'H' else ('human_lambda' if chain.chain_type == 'L' else 'human_kappa')
    html = None
    for retry in range(5):
        url = f'http://www.bioinf.org.uk/abs/shab/shab.cgi?aa_sequence={seq}&DB={chain_type}'
        request = requests.get(url)
        time.sleep(0.5 + retry * 5)
        if request.ok:
            html = request.text
            break
        else:
            print('Retry', retry+1)
    if not html:
        raise ValueError('Z-score server is not accessible')
    matches = SCORE_REGEX.findall(html)
    if not matches:
        print(html)
        raise ValueError(f'Error calling url {url}')
    return float(matches[0])

def get_z_scores_online(queries):
    results = []
    for query in tqdm(queries):
        zscore = get_z_score_online(query.seq)
        results.append(OrderedDict(
            id=query.id,
            description=query.description,
            zscore=zscore
        ))
    return pd.DataFrame(results)

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="Query FASTA sequence file.")
    parser.add_argument("output", help="Output TSV file path.")

    options = parser.parse_args()
    
    queries = SeqIO.parse(options.query, 'fasta')
    
    print('Note: The sequences will be processed through UCL Z-score web service! Sleeping for 10s, press Ctrl+C to cancel...')
    time.sleep(10)
    print('Processing...')
    df = get_z_scores_online(queries)
    
    df.to_csv(options.output, sep='\t', index=False)
    print(f'Saved to: {options.output}')
