def log2(x):
    o = 0
    while x > 1:
        x //= 2
        o += 1
    return o

def is_power_of_2(x):
    return x > 0 and x&(x-1) == 0

def raw_mul(a, b):
    if a*b == 0:
        return 0
    o = 0
    for i in range(log2(b) + 1):
        if b & (1<<i):
            o ^= a<<i
    return o

def raw_mod(a, b):
    blog = log2(b)
    alog = log2(a)
    while alog >= blog:
        if a & (1<<alog):
            a ^= (b << (alog - blog))
        alog -= 1
    return a

class BinaryField():
    def __init__(self, modulus):
        self.modulus = modulus
        self.height = log2(self.modulus)
        self.order = 2**self.height - 1
        for base in range(2, min(modulus - 1, 80)):
            powers = [1]
            while (len(powers) == 1 or powers[-1] != 1) and len(powers) < self.order + 2:
                powers.append(raw_mod(raw_mul(powers[-1], base), self.modulus))
            powers.pop()
            if len(powers) == self.order:
                self.cache = powers + powers
                self.invcache = [None] * (self.order + 1)
                for i, p in enumerate(powers):
                    self.invcache[p] = i
                return
        raise Exception("Bad modulus")

    def add(self, x, y):
        return x ^ y

    sub = add

    def mul(self, x, y):
        return 0 if x*y == 0 else self.cache[self.invcache[x] + self.invcache[y]]

    def sqr(self, x):
        return 0 if x == 0 else self.cache[(self.invcache[x] * 2) % self.order]

    def div(self, x, y):
        return 0 if x == 0 else self.cache[self.invcache[x] + self.order - self.invcache[y]]

    def inv(self, x):
        assert x != 0
        return self.cache[(self.order - self.invcache[x]) % self.order]

    def exp(self, x, p):
        return 1 if p == 0 else 0 if x == 0 else self.cache[(self.invcache[x] * p) % self.order]

    def multi_inv(self, values):
        partials = [1]
        for i in range(len(values)):
            partials.append(self.mul(partials[-1], values[i] or 1))
        inv = self.inv(partials[-1])
        outputs = [0] * len(values)
        for i in range(len(values), 0, -1):
            outputs[i-1] = self.mul(partials[i-1], inv) if values[i-1] else 0
            inv = self.mul(inv, values[i-1] or 1)
        return outputs

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    # Evaluate a polynomial at a point
    def eval_poly_at(self, p, x):
        y = 0
        power_of_x = 1
        for i, p_coeff in enumerate(p):
            y ^= self.mul(power_of_x, p_coeff)
            power_of_x = self.mul(power_of_x, x)
        return y
        
    # Arithmetic for polynomials
    def add_polys(self, a, b):
        return [((a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0))
                for i in range(max(len(a), len(b)))]

    sub_polys = add_polys
    
    def mul_by_const(self, a, c):
        return [self.mul(x, c) for x in a]
    
    def mul_polys(self, a, b):
        o = [0] * (len(a) + len(b) - 1)
        for i, aval in enumerate(a):
            for j, bval in enumerate(b):
                o[i+j] ^= self.mul(a[i], b[j])
        return o
    
    def div_polys(self, a, b):
        assert len(a) >= len(b)
        a = [x for x in a]
        o = []
        apos = len(a) - 1
        bpos = len(b) - 1
        diff = apos - bpos
        while diff >= 0:
            quot = self.div(a[apos], b[bpos])
            o.insert(0, quot)
            for i in range(bpos, -1, -1):
                a[diff+i] ^= self.mul(b[i], quot)
            apos -= 1
            diff -= 1
        return o

    # Build a polynomial that returns 0 at all specified xs
    def zpoly(self, xs):
        root = [1]
        for x in xs:
            root.insert(0, 0)
            for j in range(len(root)-1):
                root[j] ^= self.mul(root[j+1], x)
        return root

def modular_inverse(x, n):
    # Fermet's little theorem
    # -- x/y ==> (x * y^p-2) % p where p can be any prime number
    # -- let 3's inverse be a.  3*a = 1. in rational number system, a=1/3. in a finite field system where p=7, a=5 because 3*5 mod 7 = 1
    return pow(x, n - 2, n)

def fft(vals, modulus, domain):
    if len(vals) == 1:
        return vals

    L = fft(vals[::2], modulus, domain[::2])
    R = fft(vals[1::2], modulus, domain[::2])
    o = [0 for i in vals]

    for i, (x, y) in enumerate(zip(L, R)):
        y_times_root = y * domain[i]
        o[i] = (x + y_times_root) % modulus
        o[i+len(L)] = (x - y_times_root) % modulus
    return o

def inverse_fft(vals, modulus, domain):
    vals = fft(vals, modulus, domain)
    return [x * modular_inverse(len(vals), modulus) % modulus for x in [vals[0]] + vals[1:][::-1]]

# return x*y but through fft algorithm
# based on 10-digit system
def mul_with_fft(x, y, modulus, domain):
    # 1. transform x and y into an array that equals in meaning a polynomial for each number
    #    this means that this array is equaivalent to a polynomial representation.
    dl = len(domain)
    px = []
    yx = []
    mask = 10
    for i in range(dl):
        p = pow(10, i)
        cx = (x // p) % mask
        cy = (y // p) % mask
        px.append(cx)
        yx.append(cy)
    
    # 2. fft for each value, which turns a polynomial into an evaluation form 
    f1 = fft(px, modulus, domain)
    f2 = fft(yx, modulus, domain)

    # 3. multiply two evaluations (generally speaking, do some operations on two evaluations, it's not limited to multiplication)
    f3 = [(v1 * v2) % modulus for v1, v2 in zip(f1, f2)]

    # 4. transforms back the evaluation into polynomial that indicates a final number calculated
    n = inverse_fft(f3, modulus, domain)

    # 5. handling carry
    for i in range(len(n)):
        carry = n [i]// 10
        mod = n[i] % 10
        n[i] = mod
        if i < len(n) - 1:
            n[i+1] += carry
    
    # 6. back from polynomial representation to real number
    rn = 0
    for i in range(len(n)):
        rn += pow(10, i) * n[i]
    return rn


def gen_multip_subgroup_for_fft(gen, modulus):
    l = []
    l.append(pow(gen, 0, modulus))  # must be 1
    i = 1

    while True:
        r = pow(gen, i, modulus)
        if r == 1 and len(l) % 2 == 0:
            return l
        elif r == 1 and len(l) % 2 != 0:
            print("cannot find group!")
            return []
        else:
            i += 1
            l.append(r)
    
    print("cannot find group!")
    return []

if __name__ == "__main__":
    modulus = 337
    gen = 85
    multip_subgroup = gen_multip_subgroup_for_fft(gen, modulus)
    px = [3,1,4,1,5,9,2,6]
    print(multip_subgroup)

    values = fft(px, modulus, multip_subgroup)
    interp_px = inverse_fft(values, modulus, multip_subgroup)

    print(values)
    print(px)
    print(interp_px)

    print(mul_with_fft(1253, 1895, modulus, multip_subgroup))
    print(1253 * 1895)