"""Computational Synteny Block (CSB) obtainer from GECKO's frags

This script allows the user to, given a file containing frags, "merge"
them in larger structures (CSBs).

This tool accepts a txt file (header of the frags file) and a comma 
separated value file (.csv), containing the actual frags.

This script requires that `pandas`, `numpy` and `argparse` be installed 
within the Python environment you are running this script in.

Loads a frags file to a dataframe, creates the diagonal column,
sorts it by its diagonal and xStart column, divides it in n files 
and runs them. Finally, the files are concatenated and the suitable.
"""

import pandas as pd
import numpy as np
import argparse
import itertools
from common_functions import update_csb, is_suitable

def frags_to_csb (df, x_len, y_len):
    it_frags = df.iterrows()
    current_frag = next(it_frags)
    current_csb = current_frag
    dropped_frags = dict()

    df_csb = pd.DataFrame()
    try:
        while True:
            it_frags, it_candidates = itertools.tee(it_frags)
            next_candidate = next(it_candidates)
            continue_value = 1
            try:
                while continue_value > 0.1 :
                    if (is_suitable(current_frag[1], next_candidate[1])):
                        dropped_frags[next_candidate[0]] = True
                        current_csb = update_csb(current_csb, next_candidate[1], max(x_len, y_len))
                        continue_value = min(continue_value * 2, 1)
                    else:
                        continue_value = continue_value / 2
                    next_candidate = next(it_candidates)
            except StopIteration:
                pass
            finally:
                del it_candidates
            df_csb = df_csb.append(current_csb[1].to_frame().T)
            condition = True
            while condition:
                next_frag = next(it_frags)
                current_csb = next_frag
                current_frag = next_frag
                condition = next_frag[0] in dropped_frags
    except StopIteration:
        pass
    finally:
        df_csb = df_csb.append(current_csb[1].to_frame().T)
        del it_frags

    return df_csb

parser = argparse.ArgumentParser(description='asd')
parser.add_argument('header_filename', type=str, nargs=1, help='Absolute or relative path to the file containing the gecko frag file header')
parser.add_argument('frags_filename', type=str, nargs=1, help='Absolute or relative path to the file containing the gecko frag file frags (csv)')

args = parser.parse_args()
header_filename = args.header_filename[0]
frags_filename = args.frags_filename[0]

x_len = 0
y_len = 0

with open(header_filename, 'r') as header_file:
    for line in header_file:
        if 'SeqX length' in line:
            x_len = int(line.split()[3])
        if 'SeqY length' in line:
            y_len = int(line.split()[3])

df = pd.read_csv(frags_filename)

df = df.loc[df['strand(f/r)'] == 'f']
df['diagonal'] = df['xStart'] - df['yStart'] // (np.sqrt(max(x_len, y_len)))

df_sort_diag_x = df.sort_values(['diagonal', 'xStart'])

#df_sort_diag_x.to_csv('sorted_frags.csv', index = False)

df_csb = frags_to_csb(df_sort_diag_x, x_len, y_len)

df_csb.to_csv('temp_csb.csv', index = False)
filenames = [header_filename, 'temp_csb.csv']
with open('csb_final.csv', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)