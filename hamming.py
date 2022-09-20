# Linear block code (e.g., Hamming and Reed Solomon Code)
import numpy as np
from bitarray import bitarray
#from matplotlib import pyplot as plt

# codeword = message * generator
def simple_codeword_gen(msg, generator):
    return msg * generator

def test_naive_lbc():
    m = np.array([1,0,1,0])
    g1 = np.array([1,0]).reshape(-1,1)
    print(simple_codeword_gen(m, g1))

def print_graph():
    data = np.array([
        [1, 2],
        [2, 3],
        [3, 6],
    ])

# assume that it works in GF(2)
class HammingCode():
    def __init__(self):
        n = 7
        k = 4
        self.n = n
        self.k = k
        self.p = n - k  # the number of parity bits, generally 3.
        self.find_parity_matrix()

    def find_parity_matrix(self):
        mlist = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]
        P = np.zeros(shape=(4,3), dtype=int)
        I = np.identity(4, dtype=int)

        # codeword (C) = message (M) * parity_matrix (G), where G = identity_matrix (I) || parity_matrix (P), P is a matrix of (4,3)
        # P = [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3], [x4,y4,z4]]

        # In (7,4) hamming code, there are a message of 4-bits and a parity of 3-bits.
        # So, 
        for i, m in enumerate(mlist):
            a1 = (m[0] + m[1] + m[2]) % 2
            a2 = (m[0] + m[1] + m[3]) % 2
            a3 = (m[1] + m[2] + m[3]) % 2
            P[i][0] = a1
            P[i][1] = a2
            P[i][2] = a3

        G = np.concatenate([I, P], 1)
        self.G = G # for encoding
        # Generator matrix (G) for encoding = I_k | P
        
        # find a matrix for decoding
        PT = np.transpose(P)
        I = np.identity(3, dtype=int)
        H = np.concatenate([PT, I], 1)
        self.H = H
        self.HT = np.transpose(self.H)

        # parity check matrix (H) for decoding = P^T | I_n-k
        # it uses syndrome, which means, when we denote a received data as 's' and the data is not corrupted, it's guaranteed 's * H_T = 0' why does it work?
    
    def encode(self, m):
        M = np.array(m).reshape(1,-1)
        R = np.matmul(M, self.G)
        for i in range(4,7):
            R[0][i] = R[0][i] % 2
        return R
    
    def syndrome(self, c):
        S = np.matmul(c, self.HT)
        for i in range(0, 3):
            S[0][i] = S[0][i] % 2
        return S

    def decode(self, c, s):
        # syndrome decoding function
        idx = 0
        for i in range(len(self.H[0])):
            if self.H[0][i] == s[0][0] and self.H[1][i] == s[0][1] and self.H[2][i] == s[0][2]:
                c[0][i] = (c[0][i] + 1) % 2
        return c

    def contain_error(self, s):
        for i in s[0]:
            if i != 0:
                return True
        return False

def add_array(a, b):
    r = [0,0,0,0]
    for i, d in enumerate(a):
        idx = len(a) - i - 1
        s = a[idx] + b[idx]
        if s == 2:
            r[idx] = 0
            if idx > 0:
                r[idx-1] += 1
        else:
            r[idx] = r[idx] + s
            if r[idx] > 1:
                r[idx] = r[idx] % 2
                if idx > 0:
                    r[idx-1] += 1
    return r

def distance(x, y):
    dist = 0
    for i, _ in enumerate(x):
        if x[i] != y[i]:
            dist += 1
    return dist

def build_hm_table():
    hm = HammingCode()
    x = [0,0,0,0]
    one = [0,0,0,1]
    last = [1,1,1,1]
    hm_table = []

    while True:
        npx = np.array(x)
        entry = hm.encode(npx)
        hm_table.append(entry[0])
        if x == last:
            break
        x = add_array(x, one)

    # compute the minimum hamming distance
    min_dist = 9999
    for i, di in enumerate(hm_table):
        for j, dj in enumerate(hm_table):
            if i != j and distance(di, dj) < min_dist:
                min_dist = distance(di, dj)

    assert(min_dist == 3)
    return hm_table

if __name__ == "__main__":
    #build_hm_table()

    hm = HammingCode()
    m = [0, 0, 1, 1]

    c = hm.encode(m)
    c[0][1] = 1  # causing an error
    print("codeword: ", c)

    s = hm.syndrome(c)
    print("syndrome: ", s)
    print("contain error: ", hm.contain_error(s))

    d = hm.decode(c, s)
    print("after correction: ", d)