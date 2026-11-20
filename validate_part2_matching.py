# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Part 2: Lemma 2.1 (mean and difference), coupling checkpoint, input-clock test (symbolic)."""
import sys, math
import sympy as sp
from common import Checker, A_of

def main():
    c = Checker("Part 2 (electric and magnetic matching)")
    x, A, kE, kB, s, dE0, dB0, dE1, dB1 = sp.symbols('x A kappa_E kappa_B s dE0 dB0 dE1 dB1', real=True)
    y = 1 / x
    YE, YB = y - A - kE, y - A - kB
    xi = sp.sqrt(YB / YE); gm = sp.sqrt(YE * YB)
    ser_xi = sp.series(xi, x, 0, 2).removeO()
    c.check("xi_out = sqrt(YB/YE) = 1 + kappa x + O(x^2), kappa = (kE - kB)/2",
            sp.simplify(ser_xi - (1 + (kE - kB) / 2 * x)) == 0)
    ser_gm = sp.series(gm - y, x, 0, 1).removeO()
    c.check("g_out^-2 = sqrt(YE YB) = y - A - kbar + O(x), kbar = (kE + kB)/2",
            sp.simplify(ser_gm - (-A - (kE + kB) / 2)) == 0)
    c.check("coupling checkpoint: common decrement A forces kbar = 0",
            sp.solve(sp.Eq(-A - (kE + kB) / 2, -A), kE)[0] == -kB)
    z = 1 + s * x
    YEz = y / z + dE0 + dE1 * (z - 1); YBz = z * y + dB0 + dB1 * (z - 1)
    cE = sp.series(YEz - y, x, 0, 1).removeO(); cB = sp.series(YBz - y, x, 0, 1).removeO()
    c.check("input-clock test: constants -s + dE(1) and s + dB(1)",
            sp.simplify(cE - (-s + dE0)) == 0 and sp.simplify(cB - (s + dB0)) == 0)
    xiz = sp.series(sp.sqrt(YBz / YEz), x, 0, 2).removeO()
    c.check("input-clock test: xi_out = 1 + (s + kappa) x + O(x^2) with kappa = (dB - dE)/2",
            sp.simplify(xiz - (1 + (s + (dB0 - dE0) / 2) * x)) == 0)
    c.check("A = 4 b0 log 2: SU(2) 0.1287559034, SU(3) 0.1931338551",
            abs(A_of(2) - 0.1287559034) < 5e-11 and abs(A_of(3) - 0.1931338551) < 5e-11)
    return c.close()

if __name__ == "__main__":
    sys.exit(main())
