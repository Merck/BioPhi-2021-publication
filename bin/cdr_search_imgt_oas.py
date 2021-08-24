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
            yield Chain(sequence=None, aa_dict=aa_dict, name=item['original_name'], scheme='imgt', chain_type=chain_type, tail='')

def evaluate_cdr_hit(query, target, result={}, same_length=False):
    """
    Compare query and target sequence positions in CDR3 and framework regions, return number of matches and % identity
    
    If previous 'result' is provided, return (None, None) if an improvement has not been found.
    """
    if same_length:
        # Only consider pairs with same sequence length
        if len(query) != len(target):
            return None, None

        # Only consider pairs with same exact CDR positions
        if len(query.cdr3_dict) != len(target.cdr3_dict) or query.cdr3_dict.keys() != target.cdr3_dict.keys():
            return None, None
        if len(query.cdr2_dict) != len(target.cdr2_dict) or query.cdr2_dict.keys() != target.cdr2_dict.keys():
            return None, None
        if len(query.cdr1_dict) != len(target.cdr1_dict) or query.cdr1_dict.keys() != target.cdr1_dict.keys():
            return None, None
    
    best_cdr_matches = result.get('num_cdr_matches')
    best_fw_matches = result.get('num_fw_matches')
    
    # Exit if we don't improve matches in CDR regions
    num_cdr_matches = get_cdr_matches(query, target)
    if best_cdr_matches is not None and num_cdr_matches < best_cdr_matches:
        return None, None

    num_fw_matches = get_fw_matches(query, target)

    # If CDRs are indecisive, exit if we don't improve number of matches in Framework regions
    if num_cdr_matches == best_cdr_matches:
        if best_fw_matches is not None and num_fw_matches < best_fw_matches:
            return None, None
    
    return num_cdr_matches, num_fw_matches

def evaluate_framework_hit(query, target, result={}, same_length=False):
    """
    Compare query and target sequence positions in framework and CDR regions, return number of matches
    
    If previous 'result' is provided, return (None, None) if an improvement has not been found.
    """
    if same_length:
        # Only consider pairs with same sequence length
        if len(query) != len(target):
            return None, None
    
    best_fw_matches = result.get('num_fw_matches')
    best_cdr_matches = result.get('num_cdr_matches')
    
    # Exit if we don't improve matches in Framework regions
    num_fw_matches = get_fw_matches(query, target)
    if best_fw_matches is not None and num_fw_matches < best_fw_matches:
        return None, None
    
    num_cdr_matches = get_cdr_matches(query, target)
    
    # If Frameworks are indecisive, exit if we don't improve number of matches in CDR regions
    if num_fw_matches == best_fw_matches:
        if best_cdr_matches is not None and num_cdr_matches < best_cdr_matches:
            return None, None
    
    return num_cdr_matches, num_fw_matches

def get_fr1_matches(query, target):
    return sum(aa == target.fr1_dict.get(pos) for pos, aa in query.fr1_dict.items())

def get_cdr1_matches(query, target):
    return sum(aa == target.cdr1_dict.get(pos) for pos, aa in query.cdr1_dict.items())

def get_fr2_matches(query, target):
    return sum(aa == target.fr2_dict.get(pos) for pos, aa in query.fr2_dict.items())

def get_cdr2_matches(query, target):
    return sum(aa == target.cdr2_dict.get(pos) for pos, aa in query.cdr2_dict.items())

def get_fr3_matches(query, target):
    return sum(aa == target.fr3_dict.get(pos) for pos, aa in query.fr3_dict.items())

def get_cdr3_matches(query, target):
    return sum(aa == target.cdr3_dict.get(pos) for pos, aa in query.cdr3_dict.items())

def get_fr4_matches(query, target):
    return sum(aa == target.fr4_dict.get(pos) for pos, aa in query.fr4_dict.items())

def get_fw_matches(query, target):
    return get_fr1_matches(query, target) + get_fr2_matches(query, target) + get_fr3_matches(query, target) + get_fr4_matches(query, target)

def get_cdr_matches(query, target):
    return get_cdr1_matches(query, target) + get_cdr2_matches(query, target) + get_cdr3_matches(query, target)

def get_matches(query, target):
    return get_fw_matches(query, target) + get_cdr_matches(query, target)


if __name__ == "__main__":
    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs='+', help="Target OAS data-unit (gzipped) JSON file path(s).")
    parser.add_argument("--query", required=True, help="Input query sequences as ANARCI CSV (IMGT-aligned) file path.")
    parser.add_argument("--output", required=True, help="Output CSV file path.")
    parser.add_argument("--limit", type=int, help="Check only first N rows in each JSON file.")
    parser.add_argument("--debug", action='store_true', help="Print out each alignment.")
    parser.add_argument("--same-length", action='store_true', help="Only consider pairs with same sequence length.")
    parser.add_argument("--framework-first", action='store_true', help="Prioritize framework identity over CDR identity.")
    options = parser.parse_args()
    
    if options.debug and not options.limit:
        raise ValueError('Only use --debug with --limit')
    
    queries = Chain.from_anarci_csv(options.query, scheme='imgt', as_series=True)
        
    evaluate_hit = evaluate_framework_hit if options.framework_first else evaluate_cdr_hit   
    
    print(f'Searching {len(queries)} antibodies in {len(options.targets)} JSON files...')
    results = {name: {} for name in queries.index}
    for json_path in options.targets:
        for hit in iterate_oas_json(json_path, limit=options.limit):
            for query in queries:
                num_cdr_matches, num_fw_matches = evaluate_hit(query, hit, results[query.name], same_length=options.same_length) 
                if options.debug:
                    print(f'{hit.name} VS {query.name}:')
                    print(hit.align(query))
                    print('CDR:', num_cdr_matches, 'FW:', num_fw_matches)
                if num_cdr_matches is None:
                    continue
                # save improvement
                results[query.name] = {
                    'num_cdr_matches': num_cdr_matches,
                    'num_fw_matches': num_fw_matches,
                    'hit_name': hit.name,
                    'hit_seq': hit.seq
                }
    
    sorted_index = [name for name in queries.index if results[name]]
    sorted_hits = [Chain(results[name]['hit_seq'], scheme='imgt', name=name) for name in sorted_index]
    table = Chain.to_dataframe(sorted_hits)
    if not table.empty:
        table.insert(1, 'num_cdr_matches', [results[name]['num_cdr_matches'] for name in sorted_index])
        table.insert(2, 'num_fw_matches', [results[name]['num_fw_matches'] for name in sorted_index])
        table.insert(3, 'num_full_matches', table['num_cdr_matches'] + table['num_fw_matches'])
        table.insert(4, 'hit_name', [results[name]['hit_name'] for name in sorted_index])
        
    table.to_csv(options.output)
    print(f'Saved {len(table)} hits to: {options.output}')
