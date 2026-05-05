#!/usr/bin/env python3
import sys

if len(sys.argv) != 2:
    raise SystemExit("Usage: get_checkpoint_from_mcf.py <mcf>")

mcf = sys.argv[1]

with open(mcf) as f:
    for line in f:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        if line.startswith("checkpoint_file"):
            key, value = line.split("=", 1)
            print(value.strip().strip('"'))
            sys.exit(0)

raise SystemExit("checkpoint_file not found in MCF")
