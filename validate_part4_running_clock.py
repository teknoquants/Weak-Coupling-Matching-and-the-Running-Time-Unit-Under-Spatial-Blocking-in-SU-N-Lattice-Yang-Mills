# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Part 4: Lemma 4.1, Corollary 4.2 (symbolic); Theorem 4.3 and Corollary 4.4 (numerical)."""
import sys, math
import numpy as np
import sympy as sp
from common import Checker, A_of

def main():
    c = Checker("Part 4 (the running time unit)")
    x, A, kap, rho, s = sp.symbols('x A kappa rho s', real=True)
    xp = x + A * x ** 2
    lhs = 1 + (kap + rho * s) * x            # F(x, 1 + s x; h) to first order
    rhs = sp.series(1 + s * xp, x, 0, 2).removeO()
    sol = sp.solve(sp.Eq(sp.expand(lhs).coeff(x, 1), sp.expand(rhs).coeff(x, 1)), s)
    c.check("Lemma 4.1: first-order matching <=> (1 - rho) s = kappa", sp.simplify(sol[0] - kap / (1 - rho)) == 0)
    F = (1 + s * x) * (1 + kap * x)           # block commuting with time dilation: rho = 1
    slope = sp.diff(sp.diff(F, s), x).subs({x: 0})
    c.check("Corollary 4.2: a block with F(x,z) = z F(x,1) has input slope rho = 1, so matching needs kappa = 0", slope == 1)
    # Theorem 4.3: random admissible trajectories
    rng = np.random.default_rng(11)
    Aval = A_of(2); worst = 0.0; ok = True
    for _ in range(4000):
        B = rng.uniform(0, 0.5); xs = rng.uniform(0.02, 0.3)
        Ap, Am = Aval + B * xs, Aval - B * xs
        if Am <= 0 or Ap * xs >= 1: continue
        sv = rng.uniform(-1, 1); M = rng.uniform(0, 1)
        d0 = 1 - abs(sv) * xs - M * xs ** 2
        if d0 <= 0: continue
        Q = 1 / (1 - Ap * xs); Cxi = (abs(sv) * Ap * Q + M * (1 + Q * Q)) / d0
        Rs = rng.uniform(-1, 1)
        xi = lambda t: 1 + sv * t + Rs * M * t * t
        xx = rng.uniform(0.001, xs)
        d = Aval + rng.uniform(-1, 1) * B * xx
        xp = xx / (1 - d * xx)
        if xp > xs: continue
        lhs_ = abs(xi(xp) / xi(xx) - 1)
        ok &= lhs_ <= Cxi * xx * xx * (1 + 1e-12)
        worst = max(worst, lhs_ / (Cxi * xx * xx))
    c.check("Theorem 4.3: |xi(g_{k+1})/xi(g_k) - 1| <= C_xi g_k^4 on 4000 random admissible steps", ok, f"max ratio {worst:.3f}")
    sv = 0.7; xx = 0.05
    ratio = (1 + sv * xx / (1 - Aval * xx)) / (1 + sv * xx) - 1
    c.check("exact one-loop step: ratio - 1 = s A x^2 / [(1 - A x)(1 + s x)]",
            abs(ratio - sv * Aval * xx ** 2 / ((1 - Aval * xx) * (1 + sv * xx))) < 1e-15)
    ok2 = True
    for sv in (-1.5, -0.5, 0.5, 1.5):
        xs = min(0.5 / abs(sv), 0.5 / Aval) * 0.999
        for xx in np.linspace(1e-4, xs, 200):
            for d in (Aval * 0.5, Aval, Aval * 1.5):
                if d * xs > 0.5: continue
                xp = xx / (1 - d * xx)
                if xp > xs: continue
                ok2 &= abs((1 + sv * xp) / (1 + sv * xx) - 1) <= 6 * Aval * abs(sv) * xx ** 2 * (1 + 1e-12)
    c.check("admissible simplified bound C_xi = 6 A |s| under B x* <= A/2, A+ x* <= 1/2, |s| x* <= 1/2", ok2)
    # Corollary 4.4: telescoping and uniform tail
    sv = 0.4; xs = 0.15; b = 1 / xs
    for J in (50, 500, 5000):
        xsq = [1 / (b + Aval * (J - k)) for k in range(J + 1)]   # k = 0 is the UV end
        prod = np.prod([(1 + sv * xsq[k + 1]) / (1 + sv * xsq[k]) for k in range(J)])
        tel = (1 + sv * xsq[J]) / (1 + sv * xsq[0])
        if J == 5000:
            c.check("Corollary 4.4: product telescopes to xi(g_J)/xi(g_0) and tends to xi(g_J)",
                    abs(prod - tel) < 1e-12 and abs(prod - (1 + sv * xs)) < 2e-3, f"{prod:.6f} vs {1+sv*xs:.6f}")
    L = 10; Am = Aval
    tail = sum((b + Am * n) ** -2 for n in range(L, 10 ** 6)); bound = 1 / (Am * (b + Am * (L - 1)))
    c.check("Corollary 4.4: integral comparison sum_{n>=L} (b + A n)^-2 <= 1/[A (b + A (L-1))]", tail <= bound, f"{tail:.4e} <= {bound:.4e}")
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
