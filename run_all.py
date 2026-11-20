# MIT License
# Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# subject to the conditions of the MIT License (see LICENSE). THE SOFTWARE IS PROVIDED "AS IS",
# WITHOUT WARRANTY OF ANY KIND.
"""Run every W2a validation script; exit 0 only if all pass."""
import subprocess, sys, os, time
SCRIPTS = ["validate_part1_time_unit.py", "validate_part2_matching.py", "validate_part3_drift.py",
           "validate_part4_running_clock.py", "validate_part5_threshold.py", "validate_numerical_tables.py", "validate_part6_composition_threshold.py"]
here = os.path.dirname(os.path.abspath(__file__)); fails = 0
for s in SCRIPTS:
    t = time.time(); r = subprocess.run([sys.executable, os.path.join(here, s)], cwd=here)
    print(f"  {'PASS' if r.returncode == 0 else 'FAIL'}  {s}  ({time.time() - t:.0f} s)"); fails += r.returncode != 0
print("ALL CHECKS PASSED" if fails == 0 else f"{fails} SCRIPT(S) FAILED")
sys.exit(1 if fails else 0)
