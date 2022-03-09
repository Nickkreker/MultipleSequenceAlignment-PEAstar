import numpy as np

from typing import List, Dict, Tuple

def load_pam250_matrix():
    '''
    Loads PAM-250 matrix from file.
    Further we assume that the gap penalty is 8.
    '''
    pam250 = {}
    with open('pam250.dat', 'r') as f:
        first_line = f.readline().split()
        for line in f:
            l_1 = line[0]
            for idx, dist in enumerate(line[1:].split()):
                l_2 = first_line[idx]
                pam250[(l_1, l_2)] = -int(dist)
    return pam250

def solve_alignment2d(str1, str2, pam250):
    '''Calculates matrix of distances and matrix of parents in 2d.'''
    l1 = len(str1)
    l2 = len(str2)
    m = np.zeros((l1 + 1, l2 + 1))
    p = np.zeros((l1 + 1, l2 + 1), dtype=str)
    INT_MAX = 2**32-1
    d = 8
    for j in range(l2 + 1):
        m[0, j] = j * 8
        p[0, j] = 'l'
    for i in range(l1 + 1):
        m[i, 0] = i * 8
        p[i, 0] = 'u'
    p[0,0] = ' '

    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            mx = INT_MAX
            mx = m[i-1, j] + 8
            p[i, j] = 'u'
            if mx > m[i-1, j-1] + pam250[(str1[i-1], str2[j-1])]:
                mx = m[i-1, j-1] + pam250[(str1[i-1], str2[j-1])]
                p[i, j] = 'd'
            if mx > m[i, j-1] + 8:
                mx = m[i, j-1] + 8
                p[i, j] = 'l'
            m[i, j] = mx
    return (m, p)

def get_alignment(str1, str2, pam250):
    '''Calculates alignment of two strings using dynamic programming.'''
    m, p = solve_alignment2d(str1, str2, pam250)
    l1 = len(str1)
    l2 = len(str2)
    t = p[-1, -1]
    pos1 = l1
    pos2 = l2
    align1 = []
    align2 = []
    while(True):
        if t == 'u':
            align1.insert(0,str1[pos1 - 1])
            align2.insert(0, '_')
            pos1 -= 1
            t = p[pos1, pos2]
        elif t == 'l':
            align1.insert(0, '_')
            align2.insert(0, str2[pos2 - 1])
            pos2 -= 1
            t = p[pos1, pos2]
        elif t == 'd':
            align1.insert(0, str1[pos1 - 1])
            align2.insert(0, str2[pos2 - 1])
            pos2 -= 1
            pos1 -= 1
            t = p[pos1, pos2]
        elif t == ' ':
            break

    return np.array((align1, align2))

def get_distances2d(seqs, pam250):
    distances2d = {}
    for i, seq1 in enumerate(seqs):
        for j, seq2 in enumerate(seqs):
            m, _ = solve_alignment2d(seq1[::-1], seq2[::-1], pam250)
            distances2d[(i, j)] = m
    return distances2d

def calc_alignment_score(alignment: np.ndarray, score_matrix: Dict[Tuple[str, str], int], gap_penalty: int = 8) -> int:
    r, c = alignment.shape
    score = 0
    for i in range(c):
        for j in range(r): 
            for k in range(r):
                if j <= k:
                    continue
                if alignment[j][i] == '_' and alignment[k][i] == '_':
                    score += 0
                if alignment[j][i] == '_' or alignment[k][i] == '_':
                    score += gap_penalty
                else:
                    score += score_matrix[(alignment[j][i], alignment[k][i])]
    return score

def make_alignment(goal, seqs: List[str]) -> np.ndarray:
    seq_lengths = list(map(len, seqs))
    alignment_length = 0
    current = goal
    while current.parent:
        current = current.parent
        alignment_length += 1

    alignment = np.empty((len(seqs), alignment_length), dtype=str)
    current = goal
    c1 = goal.coords
    t = 0
    while current.parent:
        c2 = current.parent.coords
        for i in range(len(seqs)):
            if c1[i] != c2[i]:
                alignment[i][-t - 1] = seqs[i][seq_lengths[i] - 1]
                seq_lengths[i] -= 1
            else:
                alignment[i][-t - 1] = '_'
        current = current.parent
        c1 = current.coords
        t += 1
    return alignment