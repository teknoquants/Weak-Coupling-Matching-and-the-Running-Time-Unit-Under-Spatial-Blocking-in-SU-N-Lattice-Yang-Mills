# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Parts 4-5: Proposition 4.5, Proposition 5.1 (endpoint clock), Theorem 5.2 ingredients."""
import sys, math
import numpy as np
from common import Checker, H_aniso, spectrum, A_of

def main():
    c = Checker("Parts 4.5-5 (lifting and threshold normalization)")
    rng = np.random.default_rng(5)
    a, z = rng.uniform(0, 0.5, 2000), rng.uniform(0, 0.5, 2000)
    c.check("Proposition 4.5: 1 - (1 - alpha)(1 - zeta) <= alpha + zeta", np.all(1 - (1 - a) * (1 - z) <= a + z + 1e-15))
    gL, zL = 2.2, 1.07
    e1, _ = spectrum(H_aniso(math.sqrt(zL) * gL, gL / math.sqrt(zL)))
    e0, _ = spectrum(H_aniso(gL, gL))
    c.check("Proposition 5.1: gap of H(g_L, z_L) = z_L x isotropic gap (one plaquette)",
            abs((e1[1] - e1[0]) - zL * (e0[1] - e0[0])) < 1e-10)
    ok = True
    for cc in (0.5, 1.0, 2.0):
        for de in (0.5, 1.0, 1.5):
            D = (4 / (math.e * cc * de)) ** (4 / de)
            g = np.linspace(1e-3, 3, 3000)
            ok &= np.all(np.exp(-cc / g ** de) <= D * g ** 4 * (1 + 1e-12))
    c.check("Theorem 5.2: exp(-c/g^delta) <= D_delta g^4, D_delta = (4/(e c delta))^(4/delta)", ok)
    c.check("Theorem 5.2: sum_{j>=0} 4^-j = 4/3", abs(sum(4.0 ** -j for j in range(60)) - 4 / 3) < 1e-15)
    A = A_of(2); gs2 = 1 / 6; Am = A / 2
    xs = [1 / (1 / gs2 + Am * n) for n in range(200000)]
    c.check("Theorem 5.2: sum g_k^4 <= g*^4 + g*^2/A- along the backward running bound",
            sum(t * t for t in xs) <= gs2 ** 2 + gs2 / Am)
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
