#!/usr/bin/env python

import pandas as pd
import os
import re
import numpy as np
import argparse
from abnumber import Chain, Position
import gzip
import json

def iterate_oas_json(path, limit=None):
    if path.endswith('.json.gz'):
        gzipped = True
    elif path.endswith('.json'):
        gzipped = False
    else:
        raise ValueError(f'Expected .json or .json.gz file, got: {path}')
    with (gzip.open(path) if gzipped else open(path)) as f:
        i = 0
        for line in f:
            if limit and i > limit:
                break
            item = json.loads(line)
            if 'seq' not in item:
                # skip metadata row
                continue
            i += 1
            data = json.loads(item['data'])
            v_germline = item['v']
            if v_germline.startswith('IGHV'):
                chain_type = 'H'
            elif v_germline.startswith('IGLV'):
                chain_type = 'L'
            elif v_germline.startswith('IGKV'):
                chain_type = 'K'
            else:
                raise ValueError(f'Invalid germline "{v_germline}": {path}')
            aa_dict = {Position.from_string(pos, chain_type=chain_type, scheme='imgt'): aa \
                       for region, region_data in data.items() for pos, aa in region_data.items()}
            yield Chain(sequence=None, aa_dict=aa_dict, name=item['original_name'], scheme='imgt', 
                        chain_type=chain_type, tail='')

def evaluate_hit(query, target, result={}, same_length=False):
    """
    Compare all query and target sequence positions, return number of matches
    
    If previous 'result' is provided, return (None, None) if an improvement has not been found.
    """
    if same_length:
        # Only consider pairs with same sequence length
        if len(query) != len(target):
            return None
    
    best_matches = result.get('num_matches')
    
    num_matches = get_matches(query, target)

    # Exit if we don't improve matches
    if best_matches is not None and num_matches < best_matches:
        return None
    
    return num_matches

def get_matches(query, target):
    alignment = query.align(target)
    return len(alignment) - alignment.num_mutations()

if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs='+', help="Target OAS data-unit (gzipped) JSON file path(s).")
    parser.add_argument("--query", required=True, help="Input query sequences as ANARCI CSV (IMGT-aligned) file path.")
    parser.add_argument("--output", required=True, help="Output CSV file path.")
    parser.add_argument("--limit", type=int, help="Check only first N rows in each JSON file.")
    parser.add_argument("--debug", action='store_true', help="Print out each alignment.")
    parser.add_argument("--same-length", action='store_true', help="Only consider pairs with same sequence length.")
    options = parser.parse_args()
    
    if options.debug and not options.limit:
        raise ValueError('Only use --debug with --limit')
    
    queries = Chain.from_anarci_csv(options.query, scheme='imgt', as_series=True)
        
    print(f'Searching {len(queries)} antibodies in {len(options.targets)} JSON files...')
    results = {name: {} for name in queries.index}
    for json_path in options.targets:
        for hit in iterate_oas_json(json_path, limit=options.limit):
            for query in queries:
                num_matches = evaluate_hit(query, hit, results[query.name], same_length=options.same_length) 
                if options.debug:
                    print(f'{hit.name} VS {query.name}:')
                    print(hit.align(query))
                    print('matches:', num_matches)
                if num_matches is None:
                    continue
                # save improvement
                results[query.name] = {
                    'num_matches': num_matches,
                    'hit_name': hit.name,
                    'hit_seq': hit.seq
                }
    
    sorted_index = [name for name in queries.index if results[name]]
    sorted_hits = [Chain(results[name]['hit_seq'], scheme='imgt', name=name) for name in sorted_index]
    table = Chain.to_dataframe(sorted_hits)
    if not table.empty:
        table.insert(1, 'num_matches', [results[name]['num_matches'] for name in sorted_index])
        table.insert(2, 'hit_name', [results[name]['hit_name'] for name in sorted_index])
        
    table.to_csv(options.output)
    print(f'Saved {len(table)} hits to: {options.output}')
