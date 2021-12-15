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

def get_seq_germline_content(seq, scheme='imgt', per_position=False):
    chain = Chain(seq, scheme=scheme)
    germline_chain = chain.find_merged_human_germline().renumber(scheme)
    aligned = chain.align(germline_chain)
    if per_position:
        return [(pos, aa == bb) for pos, (aa, bb) in aligned]
    num_identical = len(aligned) - aligned.num_mutations()
    return num_identical / len(chain)

def get_seqs_germline_content(queries, scheme='imgt', per_position=False):
    results = []
    for query in tqdm(queries):
        if per_position:
            for pos, is_germline in get_seq_germline_content(query.seq, scheme=scheme, per_position=True):
                results.append({
                    'id': query.id,
                    'description': query.description,
                    scheme+'_pos': pos,
                    'is_germline': is_germline
                })
        else:
            results.append(dict(
                id=query.id,
                description=query.description,
                germline_content=get_seq_germline_content(query.seq)
            ))
    return pd.DataFrame(results)

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="Query FASTA sequence file.")
    parser.add_argument("output", help="Output TSV file path.")
    parser.add_argument("--scheme", default='imgt', help="Numbering scheme.")
    parser.add_argument("--per-position", action="store_true", help="Include one row for each position.")

    options = parser.parse_args()
    
    queries = SeqIO.parse(options.query, 'fasta')
    
    df = get_seqs_germline_content(queries, scheme=options.scheme, per_position=options.per_position)
    
    df.to_csv(options.output, sep='\t', index=False)
    print(f'Saved to: {options.output}')
