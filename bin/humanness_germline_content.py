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

def get_seq_germline_content(seq):
    chain = Chain(seq, scheme='imgt')
    germline_chain = chain.find_merged_human_germline()
    aligned = chain.align(germline_chain)
    num_identical = len(aligned) - aligned.num_mutations()
    return num_identical / len(chain)

def get_seqs_germline_content(queries):
    results = []
    for query in tqdm(queries):
        germline_content = get_seq_germline_content(query.seq)
        results.append(OrderedDict(
            id=query.id,
            description=query.description,
            germline_content=germline_content
        ))
    return pd.DataFrame(results)

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="Query FASTA sequence file.")
    parser.add_argument("output", help="Output TSV file path.")

    options = parser.parse_args()
    
    queries = SeqIO.parse(options.query, 'fasta')
    
    df = get_seqs_germline_content(queries)
    
    df.to_csv(options.output, sep='\t', index=False)
    print(f'Saved to: {options.output}')
