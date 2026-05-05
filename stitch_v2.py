#!/usr/bin/env python3
"""
stitch_all.py — corrected full‑domain‑safe version

Changes:
  - Use stitched x/y (full-domain) instead of tile00 x/y
  - Remove any xarray-added _FillValue / missing_value attributes from coordinates
  - Ensure x/y contain no NaNs or masked values
  - Fail early if coordinates look malformed
  - Preserve stitching and compression logic as before

# The content of this file has been produced with the assistance of
# University of Leeds, m365 co-pilot
"""

import argparse
import os
import sys
import subprocess
import re
import time
import shutil
from pathlib import Path
from collections import defaultdict


def msg(s: str):
    print(s, flush=True)


# ----------------------------------------------
# Utility process launcher
# ----------------------------------------------

def popen_with_log(argv, logfile: Path, extra_env=None):
    logfile.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)

    # nccopy threading
    nccopy_threads = env.get("NCCOPY_OMP_THREADS")
    if nccopy_threads:
        env["OMP_NUM_THREADS"] = str(nccopy_threads)
    else:
        env.setdefault("OMP_NUM_THREADS", str(os.cpu_count()))

    # HDF5 behaviour
    env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    env.setdefault("HDF5_PLUGIN_PATH", "")
    env.setdefault("HDF5_PARAPROCESS", "0")

    f = open(logfile, "wb")
    return subprocess.Popen(argv, stdout=f, stderr=subprocess.STDOUT, env=env)


_TILE_RE = re.compile(r"\.x(\d+)\.y(\d+)\.", re.IGNORECASE)


def decode_xy(fname):
    m = _TILE_RE.search(fname)
    if not m:
        raise ValueError(f"Cannot decode tile indices from '{fname}'")
    return int(m.group(1)), int(m.group(2))


def group_files_by_tile(files):
    out = defaultdict(list)
    for p in files:
        x, y = decode_xy(p.name)
        out[(x, y)].append(p)
    for k in out:
        out[k].sort()
    return out


# ---------------------------------------------------------
# Coordinate cleaning
# ---------------------------------------------------------

def clean_stitched_coordinates(ds_out):
    """
    Clean xarray-stitched coordinate variables.

    We keep the stitched full-domain x and y but:
      - ensure they are plain numpy arrays (no masks)
      - ensure they contain no NaNs
      - remove any _FillValue / missing_value attributes
    """

    import numpy as np

    for coord in ["x", "y"]:
        if coord not in ds_out.coords:
            continue

        vals = ds_out[coord].values
        vals = np.asarray(vals)  # drop masks

        if np.isnan(vals).any():
            raise SystemExit(
                f"[FATAL] stitched coordinate '{coord}' contains NaNs; "
                "tiles missing or stitching misalignment."
            )

        ds_out = ds_out.assign_coords({coord: (coord, vals)})

        for bad in ["_FillValue", "missing_value", "fill_value"]:
            ds_out[coord].attrs.pop(bad, None)

    return ds_out


# ---------------------------------------------------------
# Stitching for a single set
# ---------------------------------------------------------

def stitch_one_set(
    dump_dir: Path,
    set_tag: str,
    out_path: Path,
    n_threads: int,
    dask_tmp: str,
    time_to_timedelta=False,
    base_time=None,
    parallel_open=False
):
    import numpy as np
    import pandas as pd
    import xarray as xr
    import dask

    dask.config.set(scheduler="threads")
    dask.config.set(num_workers=n_threads)

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    os.environ["DASK_NUM_WORKERS"] = str(n_threads)

    if dask_tmp:
        os.environ["DASK_TEMPORARY_DIRECTORY"] = dask_tmp

    # Collect tiles
    rows = sorted([p for p in dump_dir.iterdir()
                   if p.is_dir() and re.fullmatch(r"row\d{3}", p.name)])
    all_nc = []
    for r in rows:
        sdir = r / set_tag
        if sdir.is_dir():
            all_nc.extend(sorted(sdir.glob("*.nc")))

    if not all_nc:
        raise SystemExit(f"No tiles found for {set_tag}")

    tile_files = group_files_by_tile(all_nc)
    xs = sorted({x for x, _ in tile_files})
    ys = sorted({y for _, y in tile_files})

    # Grid of files
    grid = [[tile_files.get((x, y), []) for x in xs] for y in ys]

    msg(f"[INFO] Opening tiles for {set_tag}")

    ds = xr.open_mfdataset(
        grid,
        engine="scipy",
        combine="nested",
        concat_dim=["y", "x", "time"],
        decode_times=False,
        parallel=parallel_open,
        chunks="auto"
    )

    # Which vars are spatial?
    stitched_vars = [
        v for v in ds.data_vars
        if ("x" in ds[v].dims and "y" in ds[v].dims)
    ]

    # Choose tile00 for global non-(x,y) vars
    x0, y0 = min(xs), min(ys)
    tile00 = xr.open_mfdataset(
        tile_files[(x0, y0)],
        engine="scipy",
        combine="by_coords",
        decode_times=False
    )

    global_vars = [
        v for v in tile00.data_vars
        if not ("x" in tile00[v].dims and "y" in tile00[v].dims)
    ]

    ds_out = ds[stitched_vars]
    if global_vars:
        ds_out = xr.merge([ds_out, tile00[global_vars]])

    # Time corrections
    if "time" in ds_out.dims:
        t = ds_out["time"]
        if np.issubdtype(t.dtype, np.number):
            if time_to_timedelta:
                td = pd.to_timedelta(t.values.astype(float), unit="s")
                ds_out = ds_out.assign_coords(time=("time", td))
            elif base_time:
                base = np.datetime64(pd.Timestamp(base_time))
                td = pd.to_timedelta(t.values.astype(float), unit="s")
                ds_out = ds_out.assign_coords(
                    time=("time", (base + td).astype("datetime64[ns]"))
                )
        try:
            ds_out = ds_out.sortby("time")
        except Exception:
            pass

    # Coordinate cleaning (important fix)
    ds_out = clean_stitched_coordinates(ds_out)

    msg(f"[INFO] Writing uncompressed NetCDF: {out_path}")
    encoding = {v: {"zlib": False} for v in ds_out.data_vars}
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")

    ds_out.to_netcdf(
        tmp,
        engine="h5netcdf",
        encoding=encoding,
    )

    # Atomic commit
    tmp.replace(out_path)

# ---------------------------------------------------------
# Set discovery
# ---------------------------------------------------------

def discover_unique_sets(dataset_dir):
    sets = set()
    if not dataset_dir.is_dir():
        return []
    for row in sorted(dataset_dir.iterdir()):
        if row.is_dir() and re.fullmatch(r"row\d{3}", row.name):
            for s in row.iterdir():
                if s.is_dir() and re.fullmatch(r"set\d{4}", s.name):
                    sets.add(s.name)
    return sorted(sets)

# ---------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", required=True)

    # Stage control
    ap.add_argument("--skip-stage-a", action="store_true")
    ap.add_argument("--skip-stage-b", action="store_true")

    # Stage A resources
    ap.add_argument("--concurrency-a", type=int, default=8)
    ap.add_argument("--threads-a", type=int, default=8)

    # Stage B resources
    ap.add_argument("--concurrency-b", type=int, default=8)
    ap.add_argument("--threads-b", type=int, default=8)

    # Run control
    ap.add_argument("--clean", action="store_true",
                    help="Ignore existing outputs and re-run requested stages")
    ap.add_argument("--keep-stage-a", action="store_true",
                    help="Do not delete Stage A outputs when Stage B completes")
    ap.add_argument("--only-set", default=None,
                    help="Restrict to sets: set0001,set0003 or set0000:set0007")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--parallel-open",
        action="store_true",
        help="Use open_mfdataset(parallel=True) in Stage A"
    ) 
    # Misc
    ap.add_argument("--logdir", default="logs_stitch")
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--as-timedelta", action="store_true")
    ap.add_argument("--base", default=None)

    # nccopy config (Stage B)
    ap.add_argument("--nccopy-k", default="netCDF4")
    ap.add_argument("--nccopy-d", type=int, default=1)
    ap.add_argument("--nccopy-m", default="512M")
    ap.add_argument("--nccopy-bin", default=None)
    args = ap.parse_args()

    root = Path(os.environ.get("PBS_O_WORKDIR", os.getcwd()))
    logdir = (root / args.logdir).resolve()
    logdir.mkdir(parents=True, exist_ok=True)

    datasets = [Path(d).resolve() for d in args.datasets]
    dask_tmp = args.tmpdir or f"/tmp/{os.environ.get('USER','u')}_dask"
    Path(dask_tmp).mkdir(parents=True, exist_ok=True)

    python_bin = shutil.which("python") or sys.executable
    ionice = shutil.which("ionice")

    nccopy = args.nccopy_bin or shutil.which("nccopy")
    if not nccopy:
        site_default = "/common/metoffice/nco/bin/5.0.6/nccopy"
        if Path(site_default).is_file() and os.access(site_default, os.X_OK):
            nccopy = site_default

    tasks = []
    for d in datasets:
        for s in discover_unique_sets(d):
            tasks.append((d, s))

    def parse_only_set(arg):
        if arg is None:
            return None

        out = set()
        arg = arg.strip()

        if ":" in arg:
            a, b = arg.split(":", 1)
            lo, hi = int(a[3:]), int(b[3:])
            for i in range(lo, hi + 1):
                out.add(f"set{i:04d}")
        else:
            for part in arg.split(","):
                out.add(part.strip())

        return out

    only_sets = parse_only_set(args.only_set)

    def throttle(active, limit):
        while True:
            still = []
            for item in active:
                p = item[-2]          # process is always second-to-last
                if p.poll() is None:
                    still.append(item)
            active[:] = still
            if len(active) < limit:
                return
            time.sleep(0.2)

    # Stage A
    active = []
    this_script = Path(__file__).resolve()

    for d, s in tasks:
        s=s.strip()
        # Operator intent
        if args.skip_stage_a:
            continue
        if only_sets is not None and s not in only_sets:
            continue

        out_unc = d / f"stitched_UNC_{s}.nc"
        out_fin = d / f"stitched_{s}.nc"        # Stage B output

        # Stage B dominates Stage A
        if out_fin.is_file() and not args.clean:
            if out_unc.is_file():
                try:
                    out_unc.unlink()
                except Exception:
                    pass
            continue

        # Stage A already complete
        if out_unc.is_file() and not args.clean:
            continue

        throttle(active, args.concurrency_a)

        log = logdir / f"stitch_{d.name}_{s}.log"

        argv = []
        if ionice:
            argv += [ionice, "-c2", "-n7"]
        argv += [
            python_bin,
            str(this_script),
            "--_stitch-one",
            "--dump-dir", str(d),
            "--set-tag", s,
            "--output", str(out_unc),
            "--threads", str(args.threads_a),
            "--tmpdir", dask_tmp,
        ]

        if args.as_timedelta:
            argv.append("--as-timedelta")
        if args.base:
            argv += ["--base", args.base]
        if args.dry_run:
            argv.append("--dry-run")
        if args.parallel_open:
            argv.append("--parallel-open")
        p = popen_with_log(argv, log)
        active.append((f"stitch:{d.name}/{s}", p, log))
    
    for _, p, _ in active:
        p.wait()
    active.clear()

    # -------------------------
    # Stage B: compression
    # -------------------------

    active = []
    this_script = Path(__file__).resolve()

    for d, s in tasks:
        # ---- operator intent ----
        if args.skip_stage_b:
            continue
        if only_sets is not None and s not in only_sets:
            continue

        # ---- file paths ----
        src = d / f"stitched_UNC_{s}.nc"
        dst = d / f"stitched_{s}.nc"

        # ---- file‑state logic ----
        # Final already exists → skip (unless --clean), but enforce dominance
        if dst.is_file() and not args.clean:
            if src.is_file() and not args.keep_stage_a:
                try:
                    src.unlink()
                except Exception:
                    pass
            continue

        # Stage A must exist to run Stage B
        if not src.is_file():
            continue

        # ---- execution control ----
        throttle(active, args.concurrency_b)

        log = logdir / f"compress_{d.name}_{s}.log"

        argv = []
        if ionice:
            argv += [ionice, "-c2", "-n7"]

        argv += [
            python_bin, str(this_script),
            "--_compress-one",
            "--src", str(src),
            "--dst", str(dst),
            "--threads", str(args.threads_b),
            "--nccopy-bin", nccopy,
            "--nccopy-k", args.nccopy_k,
            "--nccopy-d", str(args.nccopy_d),
            "--nccopy-m", args.nccopy_m,
        ]

        if args.clean:
            argv.append("--clean")
        if args.dry_run:
            argv.append("--dry-run")
        if args.keep_stage_a:
            argv.append("--keep-stage-a")

        p = popen_with_log(argv, log)
        active.append((f"compress:{d.name}/{s}", p, log))

    # Wait for all compression workers
    for _, p, _ in active:
        p.wait()

    Path("STITCH.ALL.DONE").touch()
    msg("ALL_DONE")

    return 0

# ---------------------------------------------------------
# Hidden subcommand
# ---------------------------------------------------------

def _stitch_one_cli():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--_stitch-one", action="store_true")
    ap.add_argument("--dump-dir", required=True)
    ap.add_argument("--set-tag", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--as-timedelta", action="store_true")
    ap.add_argument("--base", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--parallel-open", action="store_true")
    args, _ = ap.parse_known_args()

    if args.dry_run:
        msg(f"[DRY] Would stitch {args.dump_dir} {args.set_tag} -> {args.output}")
        return 0

    try:
        stitch_one_set(
            dump_dir=Path(args.dump_dir),
            set_tag=args.set_tag,
            out_path=Path(args.output),
            n_threads=args.threads,
            dask_tmp=args.tmpdir,
            time_to_timedelta=args.as_timedelta,
            base_time=args.base,
            parallel_open=args.parallel_open,
        )
        return 0
    except Exception as e:
        msg(f"ERROR: stitch failed: {e}")
        return 1

def _compress_one_cli():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--_compress-one", action="store_true")
    ap.add_argument("--src", required=True,
                    help="Stage A input file (stitched_UNC_<set>.nc)")
    ap.add_argument("--dst", required=True,
                    help="Stage B output file (stitched_<set>.nc)")
    ap.add_argument("--nccopy-bin", default="nccopy")
    ap.add_argument("--nccopy-k", default="netCDF4")
    ap.add_argument("--nccopy-d", type=int, default=1)
    ap.add_argument("--nccopy-m", default="512M")
    ap.add_argument("--threads", type=int, default=None)  # optional, if you want env control
    ap.add_argument("--clean", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--keep-stage-a", action="store_true",
                    help="Do not delete Stage A outputs after compression")
    args, _ = ap.parse_known_args()

    src = Path(args.src)
    dst = Path(args.dst)
    tmp = dst.with_suffix(dst.suffix + ".tmp")

    # Idempotence: skip only if final exists and not a clean run
    if dst.is_file() and not args.clean:
        if src.is_file() and not args.keep_stage_a:
            try:
                src.unlink()
            except Exception:
                pass
        return 0

    if args.threads:
        os.environ["NCCOPY_OMP_THREADS"] = str(args.threads)

    if args.dry_run:
        msg(f"[DRY] nccopy {src} -> {dst}")
        return 0

    # Stage A must exist
    if not src.is_file():
        raise RuntimeError(f"Stage A input missing: {src}")

    # Build nccopy command
    argv = [
        args.nccopy_bin,
        "-k", args.nccopy_k,
        "-d", str(args.nccopy_d),
        "-s",
        "-m", args.nccopy_m,
        "-c", "time/1",
        str(src),
        str(tmp),
    ]

    # Run compression (blocking)
    subprocess.check_call(argv)

    # Atomic publish
    tmp.replace(dst)


    # Cleanup Stage A unless explicitly kept
    if src.is_file() and not args.keep_stage_a:
        try:
            src.unlink()
        except Exception:
            pass

    return 0

if __name__ == "__main__":
    if "--_stitch-one" in sys.argv:
        sys.exit(_stitch_one_cli())
    if "--_compress-one" in sys.argv:
        sys.exit(_compress_one_cli())
    sys.exit(main())


"""
TODO: Pipeline refactor and optional row-wise stitching
======================================================

Context
-------
Current implementation is correct and performant enough for production.
Remaining work is about *maintainability* and *structural scalability*, not
bug fixes.

Decisions explicitly taken:
- Keep flat Stage A stitching as the default.
- Row-wise stitching deferred because complexity outweighs benefit under
  current filesystem variability.
- Dask threading is system-dependent and should remain optional.

This TODO captures *intent* so future changes are deliberate.


1) Refactor orchestration (low risk, high clarity)
--------------------------------------------------

Goal:
  Improve readability and containment of responsibilities without
  changing behavior.

Concrete steps:
  - Extract Stage A orchestration into run_stage_a(...)
  - Extract Stage B orchestration into run_stage_b(...)
  - Keep worker entry points (_stitch_one_cli, _compress_one_cli) unchanged
  - Introduce a small Context/Settings object to avoid passing many args

Non-goals:
  - Do NOT introduce classes for pipeline stages
  - Do NOT introduce a workflow framework
  - Do NOT change concurrency or cleanup semantics

Rationale:
  Current `main()` is long but linear. Refactoring should *fence* logic,
  not fragment it. This will make later algorithmic changes (e.g. row-wise
  stitching) easier to reason about.


2) Optional row-wise Stage A stitching (higher complexity)
----------------------------------------------------------

Goal:
  Reduce metadata fan-in for large grids (e.g. 32x32 tiles) by stitching
  rows independently, then combining rows.

Proposed design (not implemented yet):
  - Add flag: --row-wise
  - Reuse stitch_one_set() with an added parameter: row=<rowXXX>
  - Stage A becomes:
        for each set:
            stitch rows concurrently (limited by concurrency-a)
            combine row outputs -> stitched_UNC_<set>.nc
            cleanup row temp files

Key constraints:
  - Do NOT increase default complexity for users
  - Do NOT create long-lived intermediate files
  - Reuse existing concurrency throttle
  - row-wise must coexist with --clean and --keep-stage-a semantics

Expected benefit:
  - Modest (~1.2–1.6x) Stage A speedup in best cases
  - Reduced metadata pressure
  - More predictable scaling when filesystem is quiet

Known risks:
  - Limited or no benefit on busy shared filesystems
  - Increased code complexity and failure modes
  - Harder debugging if introduced prematurely

Decision:
  Implement ONLY after orchestration refactor, so row-wise logic lives
  inside run_stage_a(), not in main().


3) Minor cleanups to consider (optional polish)
-----------------------------------------------

- Default threads-a = 1 (threads rarely help on Met Office HPC)
- Centralize Stage A deletion logic (maybe_delete_stage_a helper)
- Add short comment block documenting Stage A/B dominance invariants


Exit criteria to revisit row-wise stitching
-------------------------------------------
Reconsider implementing row-wise stitching if:
  - Stage A consistently dominates end-to-end runtime
  - Flat stitching scaling plateaus even on quiet FS
  - Grid size increases (e.g. >32x32)
  - Filesystem behavior becomes more predictable

Until then, current code is considered "good enough".

"""
