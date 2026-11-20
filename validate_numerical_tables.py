# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Numerical inputs: the P2/P4 tables (SU(2), SU(3)), the gamma check for P4, the quartic limit."""
import sys, math
import mpmath as mp
from common import Checker, A_of

TABLE = {
 2: {10: [1.086307, 1.012675, 1.179330, 1.025495, 1.387423, 1.051579],
     100: [1.567460, 1.045927, 2.451480, 1.093930, 5.957165, 1.196539],
     1000: [3.357451, 1.064522, 11.237588, 1.133171, 124.740558, 1.283919],
     10000: [8.078744, 1.067272, 65.055712, 1.139034, 4178.384429, 1.297238]},
 3: {10: [1.082513, 1.011664, 1.171165, 1.023451, 1.368540, 1.047400],
     100: [1.457941, 1.034224, 2.122053, 1.069597, 4.473511, 1.143943],
     1000: [2.484162, 1.043458, 6.158058, 1.088781, 37.605807, 1.185346],
     10000: [4.476258, 1.044675, 19.993505, 1.091323, 396.318182, 1.190886]},
}

def prods(N, K, k):
    A = A_of(N); g2 = [1 / (6 + A * (K - j)) for j in range(K + 1)]
    return (math.exp(sum(math.log1p(k * x) for x in g2)), math.exp(sum(math.log1p(k * x * x) for x in g2)))

def main():
    c = Checker("Numerical inputs (product tables)")
    for N in (2, 3):
        for K, ref in TABLE[N].items():
            got = []
            for k in (0.05, 0.10, 0.20):
                got += list(prods(N, K, k))
            c.check(f"SU({N}) K = {K}: six table entries to 6 decimals", all(abs(round(g, 6) - r) < 1.5e-6 for g, r in zip(got, ref)))
    A = A_of(2); q = 6 / A; k = 0.10; t = math.sqrt(k) / A
    for K in (10, 1000):
        exact = mp.gamma(q) ** 2 / abs(mp.gamma(q + 1j * t)) ** 2 * abs(mp.gamma(q + K + 1 + 1j * t)) ** 2 / mp.gamma(q + K + 1) ** 2
        c.check(f"P4 gamma formula at K = {K}", abs(float(exact) / prods(2, K, k)[1] - 1) < 1e-10)
    lim = float(mp.gamma(q) ** 2 / abs(mp.gamma(q + 1j * t)) ** 2)
    c.check("quartic limit 1.139718 (SU(2), kappa = 0.10)", abs(lim - 1.139718) < 5e-7, f"{lim:.6f}")
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
