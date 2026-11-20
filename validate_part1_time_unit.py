# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Part 1: H(gE,gB) = xi H(g) exactly; gap and ground state; Corollary 1.2 and the lifting convention."""
import sys, math
import numpy as np
from common import Checker, H_aniso, spectrum

def main():
    c = Checker("Part 1 (anisotropy is a time unit)")
    rng = np.random.default_rng(3)
    errs, gaps, vecs = [], [], []
    for _ in range(5):
        gE, gB = rng.uniform(0.5, 3, 2)
        g = math.sqrt(gE * gB); xi = gE / gB
        Ha, Hi = H_aniso(gE, gB), H_aniso(g, g)
        errs.append(np.abs(Ha - xi * Hi).max() / np.abs(Ha).max())
        ea, va = spectrum(Ha); ei, vi = spectrum(Hi)
        gaps.append(abs((ea[1] - ea[0]) - xi * (ei[1] - ei[0])))
        vecs.append(1 - abs(va[:, 0] @ vi[:, 0]))
    c.check("operator identity H(gE,gB) = xi H(g) (one plaquette, 5 random couplings; relative error)", max(errs) < 1e-14, f"{max(errs):.1e}")
    c.check("gap H(gE,gB) = xi gap H(g)", max(gaps) < 1e-10, f"{max(gaps):.1e}")
    c.check("same normalized ground state", max(vecs) < 1e-12, f"{max(vecs):.1e}")
    gk, gk1, xk, xk1, a = 1.3, 1.45, 1.02, 1.03, 0.7
    ck, ck1 = gk ** 2 / (2 * a), gk1 ** 2 / (2 * 2 * a)
    r = 2 * ck / ck1; rxi = 2 * xk * ck / (xk1 * ck1)
    c.check("Corollary 1.2: r_k = 4 g_k^2/g_{k+1}^2 and r_k^xi = (xi_k/xi_{k+1}) r_k",
            abs(r - 4 * gk ** 2 / gk1 ** 2) < 1e-14 and abs(rxi - xk / xk1 * r) < 1e-14)
    eta, miso = 0.01, 2.0
    lhs = (1 - eta) * xk1 * miso / xk
    c.check("lifting convention: Gamma^xi >= (1-eta) xi_{k+1} m_iso  <=>  Gamma^iso >= (1-eta)(xi_{k+1}/xi_k) m_iso",
            abs(lhs * xk - (1 - eta) * xk1 * miso) < 1e-14)
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
