#!/usr/bin/env python3
import sys
import os

# The content of this file has been produced with the assistance of
# University of Leeds, m365 co-pilot

if len(sys.argv) != 2:
    raise SystemExit("Usage: check_continuation.py <logfile>")

logfile = sys.argv[1]

if not os.path.isfile(logfile):
    raise SystemExit(f"Logfile not found: {logfile}")

with open(logfile) as f:
    text = f.read()

if not text.strip():
    raise SystemExit(f"Logfile is empty: {logfile}")

STOP_PATTERNS = [
    "Model run complete due to model time",
    "messages file containing termination command",
    "timestep completion",
]

for pat in STOP_PATTERNS:
    if pat in text:
        print("stop")
        sys.exit(0)

print("continue")
