import numpy as np

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


def get_alignment(str1, str2, pam250, print_res=False):
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

    print(m)
    print(p)
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

    return (align1, align2)
