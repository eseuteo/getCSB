"""Common Functions for getCSB

This script contains several functions to be used by the getCSB 
script.

This file can be imported as a module and contains the following
functions:

    * is_overlapping(frag_a, frag_b) - returns true if the two
        frags overlap
    * get_distance(frag_a, frag_b) - returns the euclidean distance
        between frag_a's end and frag_b's starting point
    * is_suitable(frag_a, frag_b) - returns if the combined score
        of two frags, minus a certain penalty, is greater than 0
    * euclidean_distance(a, b) - returns the euclidean distance
        between two points a and b
    * update_csb(csb, frag, longest_seq_len) - returns a new csb 
        resulting from merging the csb and the frag. The argument
        longest_seq_len is just a way of easily obtaining the diagonal
        of the returned CSB.
"""

import numpy as np
import pandas as pd

PENALTY = .5

def is_overlapping(frag_a, frag_b):
    return frag_a['xEnd'] > frag_b['xEnd']

def get_distance(frag_a, frag_b):
    end_point = (frag_a['xEnd'], frag_a['yEnd'])
    starting_point = (frag_b['xStart'], frag_b['yStart'])
    return euclidean_distance(end_point, starting_point)

def is_suitable(frag_a, frag_b):
    if not is_overlapping(frag_a, frag_b):
        distance = get_distance(frag_a, frag_b)
        return frag_a['score'] + frag_b['score'] - distance * PENALTY > 0 
    return False

def euclidean_distance(a, b):
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)

def update_csb(csb, frag, longest_seq_len):
    distance = np.round(get_distance(csb[1], frag))
    csb[1]['xEnd'] = frag['xEnd']
    csb[1]['yEnd'] = frag['yEnd']
    csb[1]['length'] = csb[1]['length'] + distance + frag['length']
    csb[1]['score'] = csb[1]['score'] + frag['score'] - distance * PENALTY
    csb[1]['diagonal'] = csb[1]['xStart'] - csb[1]['yStart'] // (np.sqrt(longest_seq_len))
    return csb