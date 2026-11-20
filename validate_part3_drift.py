# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Part 3: Theorem 3.1 (gamma product and power law), Corollary 3.2, running checkpoint."""
import sys, math
import mpmath as mp
from common import Checker, A_of

def P2(K, b, A, k):
    return math.exp(sum(math.log1p(k / (b + A * n)) for n in range(K + 1)))

def main():
    c = Checker("Part 3 (drift of the isotropic reference)")
    A = A_of(2); b = 6.0; q = b / A
    ok = True
    for k in (0.05, 0.1, -0.1):
        s = k / A
        for K in (10, 100, 1000):
            exact = mp.gamma(K + 1 + q + s) * mp.gamma(q) / (mp.gamma(K + 1 + q) * mp.gamma(q + s))
            ok &= abs(P2(K, b, A, k) / float(exact) - 1) < 1e-10
    c.check("Theorem 3.1: product equals the gamma-function formula (kappa = 0.05, 0.1, -0.1)", ok)
    k = 0.1; s = k / A
    pref = float(mp.gamma(q) / mp.gamma(q + s))
    rat = [P2(K, b, A, k) / (pref * K ** s) for K in (10 ** 3, 10 ** 5, 10 ** 6)]
    c.check("Theorem 3.1: P2(K) ~ [Gamma(q)/Gamma(q+sigma)] K^(kappa/A)", abs(rat[-1] - 1) < abs(rat[0] - 1) and abs(rat[-1] - 1) < 2e-3,
            str([round(r, 6) for r in rat]))
    c.check("Corollary 3.2: exponent kappa/A = 0.776663 for SU(2), kappa = 0.10", abs(k / A - 0.776663) < 5e-7, f"{k/A:.6f}")
    c.check("Corollary 3.2: kappa = 0 removes the drift; kappa < 0 gives decay",
            P2(10000, b, A, 0.0) == 1.0 and P2(10000, b, A, -0.05) < P2(100, b, A, -0.05) < 1)
    # running checkpoint: y_{n+1} - y_n = A + O(1/y_n) gives y_n = b + A n + O(log(n+2))
    y = b; dev = []
    for n in range(200000):
        y = y + A + 0.3 / y
        if n in (999, 9999, 199999):
            dev.append((y - (b + A * (n + 1))) / math.log(n + 3))
    c.check("running checkpoint: y_n - (b + A n) = O(log(n+2))", max(abs(d) for d in dev) < 0.3 / A * 1.2, str([round(d, 4) for d in dev]))
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
