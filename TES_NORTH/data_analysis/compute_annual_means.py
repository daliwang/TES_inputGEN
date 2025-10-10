import argparse
import datetime as dt
import re
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import xarray as xr
import glob as _glob


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute annual means of time-series variables from an ELM NetCDF file.",
    )
    parser.add_argument(
        "--in",
        dest="input",
        required=False,
        default=None,
        help="Path to input NetCDF file (e.g., uELM_...elm.h1.0006-01-01-00000.nc). Optional if using --monthly-glob/--monthly-files",
    )
    parser.add_argument(
        "--out",
        dest="output",
        default=None,
        help="Path to output NetCDF file. If omitted, a default name is derived.",
    )
    parser.add_argument(
        "--vars",
        dest="vars",
        default=None,
        help="Comma-separated variable names to include (default: all numeric vars with time dim)",
    )
    parser.add_argument(
        "--drop-time-dim",
        action="store_true",
        help="If set, drop the time dimension in the output (single snapshot).",
    )
    parser.add_argument(
        "--chunks",
        dest="chunks",
        default="auto",
        help=(
            "Dask chunking spec. Use 'auto' (default) or a comma-separated list like "
            "'time:168,lat:100,lon:100'. Use '--no-dask' to disable."
        ),
    )
    parser.add_argument(
        "--no-dask",
        dest="no_dask",
        action="store_true",
        help="Disable dask-backed reads (may increase memory usage).",
    )
    parser.add_argument(
        "--engine",
        dest="engine",
        default=None,
        help="xarray engine to use (e.g., 'netcdf4', 'h5netcdf'). Default: let xarray choose.",
    )
    parser.add_argument(
        "--monthly-glob",
        dest="monthly_glob",
        default=None,
        help=(
            "Glob pattern for monthly h0 files for a single year (e.g., "
            "'/path/...elm.h0.0005-*.nc'). If set, annual mean is computed across these files."
        ),
    )
    parser.add_argument(
        "--monthly-files",
        dest="monthly_files",
        default=None,
        help=(
            "Comma-separated list of monthly h0 files for a single year. Useful when globbing is inconvenient."
        ),
    )
    return parser.parse_args()


def derive_default_output_path(input_path: str, year: int) -> str:
    p = Path(input_path)
    name = p.name
    # Try to match prefix like: <prefix>.hX.<year>...
    m = re.match(r"^(.*?\.h[0-9]\.)\d{4}.*?\.nc$", name)
    if m:
        prefix = m.group(1)
        return str(p.with_name(f"{prefix}{year:04d}_annual.nc"))
    # Fallback: append _annual before extension
    return str(p.with_name(p.stem + "_annual.nc"))


def select_variables(ds: xr.Dataset, explicit_vars: Optional[str]) -> List[str]:
    if explicit_vars:
        requested = [v.strip() for v in explicit_vars.split(",") if v.strip()]
        existing = [v for v in requested if v in ds.data_vars]
        missing = [v for v in requested if v not in ds.data_vars]
        if missing:
            raise ValueError(f"Requested variables not found in dataset: {', '.join(missing)}")
        return existing
    # Default: all numeric variables that have a time dimension
    var_names: List[str] = []
    for var_name, data_array in ds.data_vars.items():
        if "time" in data_array.dims and np.issubdtype(data_array.dtype, np.number):
            var_names.append(var_name)
    if not var_names:
        raise ValueError("No numeric variables with a 'time' dimension found in the dataset.")
    return var_names


def compute_annual_mean(ds: xr.Dataset, variables: List[str], keep_time_dim: bool) -> xr.Dataset:
    # Compute mean over time for selected variables
    reduced = ds[variables].mean(dim="time", skipna=True, keep_attrs=True)

    # Annotate variables with CF cell_methods
    for var in variables:
        cell_methods = reduced[var].attrs.get("cell_methods", "").strip()
        cell_methods = (cell_methods + " ").strip() if cell_methods else ""
        reduced[var].attrs["cell_methods"] = (cell_methods + "time: mean (annual)").strip()

    if keep_time_dim:
        # Keep a single-length time dimension using the first time value as representative
        # Only materialize a single scalar to avoid loading the full time coord
        t0 = ds["time"].isel(time=0).compute().values
        reduced = reduced.expand_dims({"time": [t0]})

    return reduced


def update_global_attrs(src: xr.Dataset, out: xr.Dataset, input_path: str) -> None:
    out.attrs.update(src.attrs)
    title_prefix = "Annual mean of time-series variables"
    if out.attrs.get("title"):
        out.attrs["title"] = f"{title_prefix} — {out.attrs['title']}"
    else:
        out.attrs["title"] = title_prefix
    timestamp = dt.datetime.utcnow().isoformat(timespec="seconds") + "Z"
    history_line = f"{timestamp}: annual mean computed from {Path(input_path).name}"
    if out.attrs.get("history"):
        out.attrs["history"] = f"{out.attrs['history']}\n{history_line}"
    else:
        out.attrs["history"] = history_line


def main() -> None:
    args = parse_args()

    # Configure decode_times with CFDatetimeCoder to avoid deprecation
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)

    # Determine chunks
    chunks: Optional[Dict[str, int]]
    if args.no_dask:
        chunks = None
    else:
        if args.chunks is None or str(args.chunks).lower() == "auto":
            chunks = "auto"  # type: ignore[assignment]
        else:
            # Parse spec like "time:168,lat:100,lon:100"
            chunks = {}
            for item in str(args.chunks).split(","):
                if not item.strip():
                    continue
                if ":" not in item:
                    raise ValueError(f"Invalid --chunks item: '{item}'. Expected name:size")
                name, size = item.split(":", 1)
                name = name.strip()
                size = int(size.strip())
                chunks[name] = size

    # Determine data source: monthly set (preferred if provided) or single file
    monthly_paths: List[str] = []
    if args.monthly_glob:
        monthly_paths.extend(sorted(_glob.glob(args.monthly_glob)))
    if args.monthly_files:
        monthly_paths.extend([p.strip() for p in args.monthly_files.split(",") if p.strip()])

    if monthly_paths:
        if len(monthly_paths) == 0:
            raise ValueError("No files matched --monthly-glob/--monthly-files")
        ds = xr.open_mfdataset(
            monthly_paths,
            decode_times=time_coder,
            chunks=chunks,
            engine=args.engine,
            combine="by_coords",
            parallel=True,
        )
        input_for_naming = monthly_paths[0]
    else:
        if not args.input:
            raise ValueError("Provide --in for single-file mode or --monthly-glob/--monthly-files for monthly mode.")
        ds = xr.open_dataset(
            args.input,
            decode_times=time_coder,
            chunks=chunks,  # enables dask-backed arrays
            engine=args.engine,
        )
        input_for_naming = args.input

    variables = select_variables(ds, args.vars)

    # Try to infer year from the time coordinate for output naming
    try:
        year_values = ds["time"].dt.year.values
        year = int(year_values[0]) if year_values.size > 0 else None
    except Exception:
        year = None

    out_path = args.output or (
        derive_default_output_path(input_for_naming, year)
        if year is not None
        else str(Path(input_for_naming).with_name(Path(input_for_naming).stem + "_annual.nc"))
    )

    annual_ds = compute_annual_mean(ds, variables, keep_time_dim=not args.drop_time_dim)
    update_global_attrs(ds, annual_ds, input_for_naming)

    # Apply basic compression to data variables
    encoding = {var: {"zlib": True, "complevel": 4, "shuffle": True} for var in variables}

    annual_ds.to_netcdf(out_path, format="NETCDF4", encoding=encoding)
    ds.close()
    annual_ds.close()

    print(f"Wrote annual means to: {out_path}")


if __name__ == "__main__":
    main()


