#!/usr/bin/env python3

"""
A script that checks if genes on a supplied list in csv format lie underneath known selective sweeeps in the Ag1000g
"""

import sys


import pandas as pd
import numpy as np
from collections import defaultdict

# Read in .csv file containing selection signals and associated metadata
signals = pd.read_csv("/home/dyern/lstm_scratch/sweeps/gambiae_sweeps_overlapping_genes.csv")

# Load the gene list of interest. In this example a list of genes showing allelic imbalance of expression in F1 cross progeny
sigup = pd.read_csv("/home/dyern/lstm_scratch/AI_metaanalysis/AI_in_at_least_4_final_220923.csv")
sweep = {}
nswept = {}

    # Loop through each signal, recording if any DE genes are found in one
for i, cols in signals.iterrows():

    if pd.isnull(cols['overlapping_genes']): # Skip signal if no overlapping genes
        continue
        
    sweptgenes = np.array(cols['overlapping_genes'].split(" "))

    # Get boolean array - if list of swept genes isin our DE genes
    overlap = np.isin(sweptgenes, sigup['ID'])
    #print(sweptgenes[overlap])
	
    sweep[str(cols['sweep_index'])] = sweptgenes[overlap]
    nswept[str(cols['sweep_index'])] = sweptgenes

genes = np.concatenate(list(sweep.values()))
print("test")
swept = sigup[np.isin(sigup['ID'], genes)]
#print(swept)   
for k,v in sweep.items():
    sweep[k] = ' '.join(v)

    # Build dataframe and add sweep metadata columns
sweptDE = pd.DataFrame.from_dict(sweep, orient='index', columns=['overlapping_DE_genes'])
sweptDE = sweptDE.reset_index().rename(columns={'index': 'Ag1000g_sweep'})
sweptDE['overlapping_genes'] = signals['overlapping_genes'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
sweptDE['chromosome'] = signals['new_contig'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
sweptDE['focus_pstart'] = signals['new_focus_pstart'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
sweptDE['focus_pstop'] = signals['new_focus_pstop'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)

wheresweep = defaultdict(dict)
whatsweep = defaultdict(list)

# Now loop through each swept gene to find name of sweeps it lies under
# And location of sweep
for gene in genes:

    for i, cols in sweptDE.iterrows():

        sweptgenes = np.array(cols['overlapping_DE_genes'].split(" "))

        if np.isin(sweptgenes, gene).any():
            wheresweep[gene]['chrom'] = cols['chromosome']
            wheresweep[gene]['focus_pstart'] = cols['focus_pstart']
            wheresweep[gene]['focus_pstop'] = cols['focus_pstop']

            whatsweep[gene].append(cols['Ag1000g_sweep'])

# Join name of sweeps column into a single string, so it fits in one column of data frame
for k,v in whatsweep.items():
    whatsweep[k] = ' '.join(v)
        
dfwhere = pd.DataFrame.from_dict(wheresweep, orient='index')
dfwhat = pd.DataFrame.from_dict(whatsweep, orient='index', columns=['Ag1000g_sweeps'])

df = pd.concat([dfwhat, dfwhere], axis=1)
df = df.reset_index().rename(columns={'index': 'ID'})

swept = swept.merge(df)

#write output to a tsv file (supply different name and directory if desired)
swept.to_csv("/home/dyern/lstm_scratch/sweeps/gene_list_selection_swept_AIin_at_least_4_uganda_phase3_031023.tsv", sep="\t")
