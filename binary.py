# binary field. i.e., GF(2^k)
# $python3 -m pip install galois  # install galois package first
# https://pypi.org/project/galois/

import galois
import numpy as np

# binary field
G=2
K=3
GEN = 2
GF = galois.GF(G**K)
S = set()

'''
for i in range(G**K):
    for j in range(G**K):
        x = GF([i,])
        y = GF([j,])
        z = x * y
        S.add(np.array(z)[0])
'''

CUR = 1
RES = ""
for i in range(G**K):
    x = GF([CUR])
    y = GF([GEN])
    z = x * y
    CUR = np.array(z)[0]
    RES += str(CUR)
    RES += ","
    S.add(CUR)

print(GF.properties)
print(RES)
print(len(S))