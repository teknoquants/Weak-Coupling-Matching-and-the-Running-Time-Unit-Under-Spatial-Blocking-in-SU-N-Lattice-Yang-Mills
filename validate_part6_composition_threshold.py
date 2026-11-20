# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Section 7: time-dilation covariance of slice-wise blocking, the physical clock and the composition
test (4.1), and the arithmetic of the threshold formalism."""
import sys, math
import numpy as np
import sympy as sp
from scipy.linalg import expm, logm, eigh
from scipy.special import ive
from common import Checker, A_of

def main():
    c = Checker("Section 7 (composition test and threshold formalism)")
    rng = np.random.default_rng(17)
    n, r = 12, 4
    M = rng.normal(size=(n, n)); H = M @ M.T + n * np.eye(n)
    P = np.zeros((n, n)); P[:r, :r] = np.eye(r); Q = np.eye(n) - P
    def feshbach(Hm, E):
        HQQ = Hm[r:, r:]; HPQ = Hm[:r, r:]
        return Hm[:r, :r] - HPQ @ np.linalg.solve(HQQ - E * np.eye(n - r), HPQ.T)
    z, E = 1.37, 0.4
    c.check("Theorem 13 (Schur/Feshbach extraction): H_eff[zH](zE) = z H_eff[H](E)",
            np.abs(feshbach(z * H, z * E) - z * feshbach(H, E)).max() < 1e-10)
    def transfer_gen(Hm, t):
        K = expm(-t * Hm)[:r, :r]
        return -logm(K).real / t
    t = 0.3
    c.check("Theorem 13 (transfer extraction): K_eff[zH](t) = z K_eff[H](z t)",
            np.abs(transfer_gen(z * H, t) - z * transfer_gen(H, z * t)).max() < 1e-9)
    x, A, c1, c2, zz = sp.symbols('x A c1 c2 z', positive=True)
    cref = lambda u: 1 + c1 * u + c2 * u ** 2
    xp = x / (1 - A * x)                                    # exact one-loop step x -> x'
    F = lambda u, zv: zv * cref(u) / cref(u / (1 - A * u))  # physical clock of the blocked theory
    xi = lambda u: 1 / cref(u)
    c.check("Theorem 15: xi(g) = 1/c_ref(g) satisfies the composition test (4.1) exactly (e_k = 0)",
            sp.simplify(F(x, xi(x)) - xi(xp)) == 0)
    ser = sp.series(F(x, 1), x, 0, 3).removeO()
    c.check("Theorem 15: the block's output clock F(x,1) = 1 + O(x^2), i.e. kappa = 0",
            sp.simplify(ser.coeff(x, 1)) == 0 and sp.simplify(ser.coeff(x, 2) + c1 * A) == 0)
    c.check("Corollary 16: the running clock has s = -c1 (first-order coefficient of 1/c_ref)",
            sp.simplify(sp.series(xi(x), x, 0, 2).removeO().coeff(x, 1) + c1) == 0)
    Aval = A_of(2); bD = math.atanh(1 / 18) / 4
    aF = ive(2, 4 * bD) / ive(1, 4 * bD); g2m = -(2 / 1.5) * math.log(aF)
    c.check("Section 7.4: beta_D = 0.013903 and g^2_match(beta_D) = 5.7010 (SU(2))",
            abs(bD - 0.0139032) < 1e-6 and abs(g2m - 5.70102) < 5e-5, f"{g2m:.5f}")
    M_minus = math.ceil((6 - 1 / g2m) / (Aval / 2)); M_one = math.ceil((6 - 1 / g2m) / Aval)
    c.check("Proposition 18: step counts M <= 91 (A- = A/2) and 46 (A- = A), SU(2), g*^2 = 1/6",
            M_minus == 91 and M_one == 46, f"{M_minus}, {M_one}")
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
