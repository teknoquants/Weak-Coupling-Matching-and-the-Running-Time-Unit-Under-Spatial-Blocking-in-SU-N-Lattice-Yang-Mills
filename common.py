# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Shared helpers for the W2a validation suite."""
import math
import numpy as np
from scipy.linalg import eigh

def A_of(N):
    """One-loop decrement of g^-2 per doubling, B2 convention: A = 4 b0 log 2, b0 = 11N/(48 pi^2)."""
    return 4 * (11 * N / (48 * math.pi ** 2)) * math.log(2)

class Checker:
    def __init__(self, name):
        self.name, self.n, self.ok = name, 0, 0
        print(f"== {name}")
    def check(self, label, cond, info=""):
        self.n += 1; self.ok += bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"  ({info})" if info else ""))
    def close(self):
        print(f"  -> {self.name}: {self.ok}/{self.n} checks passed")
        return 0 if self.ok == self.n else 1

def plaquette_TV(nmax=40):
    """SU(2) one plaquette, class functions of W in the character basis chi_n.
    T = sum over the 4 links of (-Delta_e): eigenvalue 4 n(n+2)/2 = 2 n(n+2).  V = N - ReTr W = 2 - X."""
    n = np.arange(nmax + 1)
    X = np.zeros((nmax + 1, nmax + 1))
    for i in range(nmax):
        X[i, i + 1] = X[i + 1, i] = 1.0
    T = np.diag(2.0 * n * (n + 2))
    V = 2 * np.eye(nmax + 1) - X
    return T, V

def H_aniso(gE, gB, nmax=40):
    T, V = plaquette_TV(nmax)
    return 0.5 * gE ** 2 * T + V / gB ** 2

def spectrum(H):
    e, v = eigh(H)
    return e, v
