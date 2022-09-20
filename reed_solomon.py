import numpy as np
import galois
from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial

def unsigned32(n):
  return n & 0xffffffffffffffff

# test: RS(n,k) with s-bit symbol
#       s=8, k=4, n=8, 2t=4, t=2 --> it can correct up to 2 symbol errors in the code word. (2byte errors)

# TODO: perhaps, we need to use GF(256) for this case?
# try with GF(256) to fall outcomes down to under 256.
class ReedSolomon():
    def __init__(self):
        self.s = 8
        self.k = 4
        self.n = 8
        self.t = int((self.n - self.k) / 2)
        self.GF = galois.GF(2**8)
        self.gen = 2
        self.G = self.find_generator()

    def find_generator(self):
        g = galois.Poly([1], field=self.GF)
        for i in range(self.n - self.k):
            t = pow(self.gen, i)
            c = galois.Poly([1, t], field=self.GF)
            g = g * c
        return g

    # encoding using lagrange interpolation
    def encode(self, m):
        assert(len(m) == self.k)

        # 1. making coordinates using the given message
        xl = [i for i, _ in enumerate(m)]
        yl = [i for i in m]
        x = self.GF(xl)
        y = self.GF(yl)

        # 2. build message polynomial
        MP = galois.lagrange_poly(x, y)

        # 3. get parity polynomial
        MPP = MP * pow((galois.Poly([1, 0], field=self.GF)), self.t * 2)  # make a room for t symbols
        PARITY = MPP % self.G  # find the remainder
        C = MPP - PARITY  # reflect the remainder

        # important note:
        # -- (1) PARITY = MPP % self.G
        # -- (2) MPP = a * self.G + PARITY (for some a)
        # -- (3) MPP - PARITY = a * self.G
        # -- (4) MPP - PARITY = C = a * self.G
        # :: with (4), we can argue that C is divisible by self.G unless C is corrupted.
        # :: for the roots of self.G (in here, 1,2,4,....) evaluation results will be 0 as self.G is down to 0. 
        # :: that is to say, evaluating C at the roots of self.G must be 0 if C is preserved.

        c = [int(i) for i in C.coeffs]
        return c

    def decode(self, c): # return True if no error by syndrome polynomial evaluations
        # 1. syndrome decoding
        # evaluating C at the roots of self.G
        C = galois.Poly(c, field=self.GF)
        for i in range(self.n - self.k):
            t = pow(self.gen, i)
            if C(t) != 0:
                return False
            # [Q] is it necessary to go over all the roots? --> i.e., 't' evaluations are needed? correct..
        return True

if __name__ == "__main__":
    rs = ReedSolomon()
    m = [0x19, 0x27, 0x3a, 0x40]
    c = rs.encode(m)

    # corruption
    c[0] += 1

    print(rs.decode(c))