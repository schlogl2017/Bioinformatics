#!/usr/bin/env python
# coding: utf-8


from collections import Counter
from Bio import SeqIO
from collections import defaultdict
import itertools
import time
import re
import matplotlib.pyplot as plt



def pattern_count(sequence, pattern):
    '''Return the count the input pattern found in to a give string.'''
    count = 0
    for i in range(0, len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            count = count + 1
    return count


def pattern_counter(sequence, pattern):
    '''Return the count the input pattern found in to a give string.'''
    return Counter([sequence[i:i+len(pattern)] for i in range(len(sequence)-len(pattern) + 1) 
                    if sequence[i:i+len(pattern)] == pattern])



def get_pattern_count_rgx(sequence, pattern):
    '''Return the count the input pattern found in to a give string.'''
    return len(re.findall(r'(?='+pattern+')', sequence))

def frequent_patterns(sequence, k):
    frequent_pats = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        pattern = sequence[i:i+k]
        frequent_pats[pattern] = frequent_pats.get(pattern, 0) + 1
        maxval = max(frequent_pats.values())
    return [k for k, v in frequent_pats.items() if v == maxval], maxval, frequent_pats



def most_frequent_patterns(sequence, k, n=1):
    """Returns the n most frequent patterns of lenght k from a input
    sequence."""
    return Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)]).most_common(n)

def reverse_complement_str(sequence):
    reverse = str.maketrans("ACGT","TGCA")
    return sequence.translate(reverse)[::-1]


def pattern_locations(sequence, pattern):
    matches = []
    n = len(sequence)
    m = len(pattern)
    for i in range(n - m + 1):
        if pattern == sequence[i:i+m]:
            matches.append(i)
    return matches


def get_pattern_starpos_rgx(sequence, pattern):
    """ Return a list of start positions of the matched pattern in the string"""
    return [match.start() for match in re.finditer(r'(?='+pattern+')', sequence)]


def pattern_positions(sequence, pattern):
    """ returns the position of all k-mers in sequence as a dictionary
    Autor: Leonard McMillan """
    k_positions = {}
    for i in range(len(sequence) - len(pattern) + 1):
        kmer = sequence[i:i+len(pattern)]
        k_positions[kmer] = k_positions.get(kmer,[])+[i]
    # combine kmers with their reverse complements
    k_pair_position = {}
    for kmer, pos in k_positions.items():
        k_rev = reverse_complement_str(kmer)
        if (kmer == k_rev):
            k_pair_position[kmer] = pos
        elif (kmer <= k_rev):
            k_pair_position[kmer] = sorted(pos + k_positions.get(k_rev, []))
        elif (k_rev < kmer):
            k_pair_position[k_rev] = sorted(k_positions.get(k_rev, []) + pos)
    return k_pair_position



def get_pattern_positions_rgx_yield(sequence, pattern):
    """ Yield a list of start positions of the matched pattern in the string"""
    positions = re.finditer(r'('+pattern+')', sequence)
    for position in positions:
        pos_lst =  position.span()[0]
        yield pos_lst


def get_pattern_positions_rgx_generator(sequence, pattern):
    """ Return a list of start positions of the matched pattern in the string"""
    return (pos.span()[0] for pos in re.finditer(r'('+pattern+')', sequence))


def pattern_to_numbers(pattern):
    """Returns a integer number representing the inputed 
    pattern"""
    k = len(pattern)
    base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    number= 0
    for i in range(k):
        number += base_index[pattern[i]]*4**(k-i-1)
    return number


def kmers(seq, k):
    """Returns a k_mer iterator"""
    return (seq[i:i+k] for i in range(len(seq)-k+1))


def number_to_pattern(number, k):
    """Returns a pattern representing the inputed 
    number"""
    nts = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    if k == 1:
        return nts[number]
    prefix_idx = number // 4
    reminder = number % 4
    prefix_pattern = number_to_pattern(prefix_idx, k - 1)
    return prefix_pattern + nts[reminder]



def compute_frequency(sequence, k):
    """Returns the number of times that each k-mer Pattern has already 
    appeared in the sequence."""
    freq_array = [0] * (4**k)
    for i, val in enumerate(sequence[:len(sequence) - (k - 1)]):
        freq_array[pattern_to_numbers(sequence[i:i + k])] += 1
    return [v for v in freq_array]



def faster_frequent_patterns(sequence, k):
    """Returns the most frequent k-mers in a string"""
    freq_pats = set()
    freq_array = compute_frequency(sequence, k)
    max_count = max(freq_array)
    for i in range(4**k -1):
        if freq_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            freq_pats.add(pattern)
    return freq_pats


def get_frequent_patterns_by_sorting(sequence, k):
    """Returns the most frequent k-mers in a string"""
    freq_pats = set()
    idx = [0] * (len(sequence) - k + 1)
    count = [0] * (len(sequence) - k + 1)
    for i in range(len(sequence) - k + 1):
        pattern = sequence[i:i+k]
        idx[i] = pattern_to_numbers(pattern)
        count[i] = 1
    sorted_idx = sorted(idx)
    for i in range(len(sequence) - k + 1):
        if sorted_idx[i] == sorted_idx[i - 1]:
            count[i] = count[i - 1] + 1
    max_count = max(count)
    for i in range(len(sequence) - k + 1):
        if count[i] == max_count:
            pattern = number_to_pattern(sorted_idx[i], k)
            freq_pats.add(pattern)
    return freq_pats


def clump_finder(sequence, k, window, times):
    """Returns the number of the times the clumps of patterns of k length appears in a sliding window
    inside a inputed sequence."""
    clumps = []
    for i in range(len(sequence) - window + 1):
        pattern = sequence[i:i+window]
        counts = {}
        for i in range(len(pattern) - k + 1):
            if pattern[i:i+k] not in counts:
                counts[pattern[i:i+k]] = 0
            counts[pattern[i:i+k]] += 1
        for mer in counts:
            if counts[mer] >= times and mer not in clumps:
                clumps.append(mer)
    return clumps



def better_clump_finder(sequence, k, window_size, num_timest):
    """Returns the number of the times the clumps of patterns of k 
    length appears in a sliding window
    inside a inputed sequence."""
    freq_patterns = []
    clumps = []
    for i in range(4 ** k):
        clumps.append(0)
    window = sequence[0:window_size]
    freq = compute_frequency(window, k)
    for i in range(4 ** k):
        if freq[i] >= num_timestt:
            clumps[i] = 1
    for i in range(1, len(sequence) - window_size):
        first_pat = sequence[i-1:i-1+k]
        index = pattern_to_numbers(first_pat)
        freq[index] = freq[index] - 1
        last_pat = sequence[i+window_siz-k:i+window_size]
        index = pattern_to_numbers(last_pat)
        freq[index] = freq[index] + 1
        if freq[index] >= num_timest:
            clumps[index] = 1
    for i in range(4**k):
        if clumps[i] == 1:
            pat = number_to_pattern(i, k)
            freq_patterns.append(pat)
    return freq_patterns


def get_base_compostion_stats(sequence, start):
    """Returns basics statistics from a sequence in itś forward and
    revese strands. Autor: Leonard McMillan"""
    half_genome = len(sequence)//2
    ter_C = start + half_genome
    if (ter_C > len(sequence)): # handle genome's circular nature
        ter_C = ter_C - len(sequence) + 1
    stats = {}
    for base in "ACGT":
        total = sequence.count(base)
        if (ter_C > start):                                   # case 1: ----start========ter---->
            forward_count = sequence[start:ter_C].count(base)
            reverse_count = total - forward_count
        else:                                                # case 2: ====ter--------start====>
            reverse_count = sequence[ter_C:start].count(base)
            forward_count = total - reverse_count
        stats[base] = (total, forward_count, reverse_count)
    return stats


def get_sequence_skew(sequence):
    """Returns the difference between the total number of 
    occurrences of G and the total number of occurrences of C in 
    the first i elements of the sequence. """
    skew = [0]
    for idx, element in enumerate(sequence):
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        elif sequence[idx] == 'C':
            skew.append(skew[idx] -1)
        else:
            skew.append(skew[idx])
    return skew


def get_minimum_skew(sequence):
    """Returns a position in a sequence minimizing the skew."""
    min_skew = []
    skew = get_sequence_skew(sequence)
    m_skew = min(skew)
    for idx in range(len(sequence) + 1):
        if skew[idx] == m_skew:
            min_skew.append(idx)
    return min_skew


def hamming_distance(sequence1, sequence2):
    """Return the HD form two inputed sequences"""
    return len([(x,y) for x,y in zip(sequence1, sequence2) if x != y])


def pattern_with_mismatches_finder(sequence, pattern, distance):
    """Return list of positions of the pattern with hamming distance d from one input sequence."""
    positions = []
    count = 0
    for i in range(len(sequence) -  len(pattern) + 1):
        pat = sequence[i:i+len(pattern)]
        if hamming_distance(pattern, pat)<= distance:
            count += 1
            positions.append(i)
    return positions



def hamming_pattern_count(sequence, pattern, distance):
    """Return the number of times the pattern with hamming distance d is found 
    in the inputed sequence."""
    return len([i for i in range(len(sequence) - len(pattern) + 1) if 
                hamming_distance(sequence[i:i+len(pattern)], pattern) <= distance])


def hamming_pattern_positions(sequence, pattern, distance):
    """Return list of positions of the pattern with hamming distance d from one input sequence."""
    return [i for i in range(len(sequence) - len(pattern) + 1) if 
                hamming_distance(sequence[i:i+len(pattern)], pattern) <= distance]



def get_neighbors(pattern, d):
    """Return list of all offsets of patterns with hamming distance d of `pattern"""
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return set(['A', 'C', 'T', 'G'])
    neighborhood = set()
    neighbors = get_neighbors(pattern[1:], d)
    for kmer in neighbors:
        if hamming_distance(pattern[1:], kmer) < d:
            for char in ['A', 'C', 'T', 'G']:
                neighborhood.add(char + kmer)
        else:
            neighborhood.add(pattern[0] + kmer)
    return sorted(list(neighborhood))


def frequent_patterns_with_d_mismatches(sequence,k,distance):
   """Return list of most frequent offsets of patterns of k 
   length and with hamming distance d from one input sequence."""
    counts = {}
    for i in range(len(sequence) - k + 1):
        for kmer in get_neighbors(sequence[i:i+k], distance):
            counts[kmer] = counts.get(kmer, 0) + 1
    max_count = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == max_count]


def frequent_patterns_with_mismatches_and_rev_complements(sequence, k, distance):
    """Return list of most frequent offsets of patterns and it's reverse
    complements of k length and with hamming distance d from one input sequence."""
    counts = {}
    for i in range(len(sequence)-k+1):
        for subsequence in [sequence[i:i+k], reverse_complement_str(sequence[i:i+k])]:
            for kmer in get_neighbors(subsequence, distance):
                counts[kmer] = counts.get(kmer, 0) + 1
    max_count = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == max_count]

