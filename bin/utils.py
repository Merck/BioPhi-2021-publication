from Bio import SeqIO
import os
import pandas as pd
import itertools
import numpy as np
import re
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib import pyplot as plt
from IPython.display import display, HTML

def heading(text, level=2):
    display(HTML(f'<h{level}>{text}</h{level}>'))

def iterate_fasta_index(records, idx):
    for i, record in enumerate(records):
        if i == idx[0]:
            yield record
            idx = idx[1:]
            if not len(idx):
                break
        
def iterate_single_fasta(path, limit=None, random=False, max_length=None, random_sample_from_limit=None):
    records = SeqIO.parse(path, 'fasta')
    if limit == 0:
        return []
    if random:
        if limit is None:
            raise ValueError('Random can only be used together with limit')
        if random_sample_from_limit:
            size = random_sample_from_limit
        else:
            # size = ...
            raise NotImplementedError('TODO: Count total number of sequences')
        limit = min(size, limit)
        random_idx = np.array(sorted(np.random.choice(size, limit, replace=False)))
        records = iterate_fasta_index(records, random_idx)
    elif limit:
        records = itertools.islice(records, limit)
        
    if max_length:
        return filter_max_length(records, max_length)

    return records
    
    
def iterate_fasta(paths, limit=None, max_length=None, random=False, random_sample_from_limit=None):
    """
    Iterate through fasta sequence file(s), return generator of records
    """
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        for record in iterate_single_fasta(path, limit=limit, max_length=max_length, random=random, random_sample_from_limit=random_sample_from_limit):
            yield record

            
def barplot(series, ax=None, limit=None, limit_agg_func='sum', remaining_label='(remaining {})', title=None, neglog=False, fmt='{}', na_label='N/A', xlogtickstep=2, padding=1.3, **kwargs):
    series = series[::-1]
    if 'color' in kwargs and isinstance(kwargs['color'], list):
        kwargs['color'] = kwargs['color'][::-1]
    if limit:
        series = series.sort_values()
        if len(series) > limit:
            remaining = series[:-limit]
            remaining_sum = pd.Series([remaining.agg(limit_agg_func)], [remaining_label.format(len(remaining))])
            top = series[-limit:]
            series = pd.concat([remaining_sum, top])
            series.name = '{} (Top {})'.format(top.name, limit)
    
    title = series.name if title is None else title
    if neglog:
        neglog_series = -np.log10(series)
        ax = neglog_series.plot.barh(ax=ax, title=title, **kwargs)
        xmax = int(np.ceil(neglog_series.max() * padding))
        ax.set_xlim([0, xmax])
        ax.set_xticks(list(range(0, xmax, xlogtickstep)))
        ax.set_xticklabels(['$\mathregular{10^{'+str(-x)+'}}$' for x in ax.get_xticks()]);
        for i, v in enumerate(series.values):
            text = fmt(v) if callable(fmt) else fmt.format(v)
            ax.text(-np.log10(v), i, ' {} '.format(text), va='center', ha='left')

    else:
        ax = series.plot.barh(ax=ax, title=title, **kwargs)
        cmin = series.min()
        cmax = series.max()
        ax.set_xlim(0 if cmin > 0 else cmin * padding, 0 if cmax < 0 else cmax * padding)
        for i, v in enumerate(series.values):
            text = fmt(v) if callable(fmt) else fmt.format(v)
            if pd.isnull(v): 
                if na_label:
                    text = na_label
                    v = 0
            ax.text(v, i, ' {} '.format(text), va='center', ha='right' if v < 0 else 'left')
    return ax

def plot_logo(counts, figsize=(30, 3), title='', color_scheme='chemistry'):
    """
    Plot motif logo using counts matrix (AAs in columns, counts in rows)
    """
    fig, ax = plt.subplots(figsize=figsize)
    counts = counts.drop('-', axis='columns')
    logo = logomaker.Logo(counts.reset_index(drop=True), width=0.95, color_scheme=color_scheme, ax=ax)
    plt.xticks(list(range(len(counts))), counts.index, rotation=90);
    ax.set_title(title)
    
def series2records(series, description=''):
    if isinstance(description, str):
        description = [description] * len(series)
    return [SeqRecord(Seq(seq), id=name, description=d) for (name, seq), d in zip(series.items(), description)]

def parse_netMHCIIpan_table(paths, rank_strong, rank_weak):
    if isinstance(paths, str):
        paths = [paths]
    assert rank_strong >= 1 or rank_weak >= 1, 'Ranks should be from 1-100'
    results = []
    for path in paths:
        result = pd.read_csv(path, sep='\t', skiprows=1, engine='python')
        with open(path, 'rt') as f:
            alleles = [allele for allele in f.readline().strip().split('\t') if allele]
            
        columns = ['Pos', 'Peptide', 'ID']
        for allele in alleles:
            columns += [f'{allele}.1-log50k', f'{allele}.nM', f'{allele}.Rank']
        columns += ['Ave', 'NB']
        result.columns = columns
        ranks = result[[f'{allele}.Rank' for allele in alleles]]
        result.insert(0, 'SB_Alleles', (ranks <= rank_strong).sum(axis=1))
        result.insert(1, 'WB+SB_Alleles', (ranks <= rank_weak).sum(axis=1))
        
        result.drop(columns=['Pos', 'NB'], inplace=True)
        results.append(result)
        
    return pd.concat(results, ignore_index=True).set_index('Peptide')