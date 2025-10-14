import argparse
import glob
import os
import subprocess
from typing import List

import netCDF4 as nc


EXCLUDED_VARS = {
    "time",
    "time_bnds",
    "mcsec",
    "date",
    "datesec",
    "x",
    "y",
    "gridID",
}


def detect_time_vars(nc_path: str) -> List[str]:
    ds = nc.Dataset(nc_path)
    try:
        names: List[str] = []
        for name, var in ds.variables.items():
            if name in EXCLUDED_VARS:
                continue
            if hasattr(var, "dimensions") and ("time" in var.dimensions):
                names.append(name)
        return names
    finally:
        ds.close()


def main() -> int:
    p = argparse.ArgumentParser(description="Batch plot all time-varying vars from annual_mean NetCDFs.")
    p.add_argument("--run-dir", default=".", help="Directory containing *.annual_mean.nc files. Default: current dir")
    p.add_argument("--mask", required=True, help="Path to mask NetCDF (with x, y, gridID).")
    p.add_argument("--out", default="./plots", help="Directory to write plots. Default: ./plots")
    p.add_argument("--time-index", type=int, default=0, help="Time index to plot (annual means usually have 1). Default: 0")
    p.add_argument("--cmap", default="viridis", help="Matplotlib colormap. Default: viridis")
    p.add_argument("--dpi", type=int, default=150, help="Output DPI. Default: 150")
    p.add_argument("--year-digits", type=int, choices=[4, 5], default=4, help="Digits for year in filename. Default: 4")
    p.add_argument(
        "--vars-to-plot",
        nargs="+",
        default=["GPP", "FSDS", "FLDS", "TBOT", "RH2M"],
        help="Space-separated variable names to plot. Default: GPP FSDS FLDS TBOT RH2m",
    )
    p.add_argument(
        "--show2d-script",
        default="/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_AOI_scripts/Show2DVariablesBatch.py",
        help="Path to Show2DVariablesBatch.py",
    )
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)

    pattern = os.path.join(args.run_dir, "*.annual_mean.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[warn] No files matched: {pattern}")
        return 0

    for nc_file in files:
        vars_to_plot = args.vars_to_plot

        prefix = os.path.splitext(os.path.basename(nc_file))[0]

        cmd = [
            "python3",
            args.show2d_script,
            "--mask",
            args.mask,
            "--data",
            nc_file,
            "--vars",
            *vars_to_plot,
            "--time-index",
            str(args.time_index),
            "--output-dir",
            args.out,
            "--prefix",
            prefix,
            "--cmap",
            args.cmap,
            "--dpi",
            str(args.dpi),
            "--year-digits",
            str(args.year_digits),
        ]

        print(f"[run] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())