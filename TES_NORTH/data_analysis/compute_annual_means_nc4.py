#!/usr/bin/env python3
"""
Compute annual means from monthly NetCDF files using netCDF4 only (no xarray),
with no-leap month-length weighting. Streams per year to minimize memory.

Assumptions:
- Each input file represents one month (time dimension length = 1)
- Data are effectively 1D or pseudo-1D (e.g., dims like (time, i, j) with j=1)

Usage example:
  compute_annual_means_nc4.py \
    \
    "/path/case.elm.h0.????-??.nc" \
    --vars TBOT FSDS GPP FLDS RH2M \
    --output "/path/annual/" -c 5
"""

import argparse
import os
import re
import sys
from glob import glob
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from netCDF4 import Dataset


NOLEAP_DAYS_IN_MONTH = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float64)


def expand_input_patterns(patterns: Sequence[str]) -> List[str]:
    files: List[str] = []
    for pat in patterns:
        expanded = sorted(glob(pat))
        if not expanded:
            if os.path.isfile(pat):
                expanded = [pat]
        files.extend(expanded)
    # Deduplicate while preserving order
    seen = set()
    unique_files: List[str] = []
    for f in files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)
    return unique_files


def group_files_by_year(files: Sequence[str]) -> Dict[int, List[str]]:
    year_to_files: Dict[int, List[str]] = {}
    for f in files:
        m = re.search(r"(\d{4})-(\d{2})\.nc$", os.path.basename(f))
        if not m:
            print(f"[warn] Could not parse year-month from filename: {f}", file=sys.stderr)
            continue
        year = int(m.group(1))
        year_to_files.setdefault(year, []).append(f)
    for y in list(year_to_files.keys()):
        year_to_files[y] = sorted(year_to_files[y])
    return year_to_files


def infer_case_stub(sample_file: Optional[str]) -> Optional[str]:
    if not sample_file:
        return None
    base = os.path.basename(sample_file)
    parts = base.split('.elm.h0.')
    if len(parts) >= 2:
        return parts[0]
    return None


def get_reference_metadata(nc_path: str) -> Tuple[Dict[str, Tuple[Tuple[str, ...], np.dtype, Optional[float], Dict[str, str]]], Dict[str, int]]:
    """
    Inspect a reference file and return:
    - var_meta: mapping var_name -> (dims, dtype, fill_value, attrs)
    - dim_sizes: mapping dim_name -> length
    """
    var_meta: Dict[str, Tuple[Tuple[str, ...], np.dtype, Optional[float], Dict[str, str]]] = {}
    dim_sizes: Dict[str, int] = {}
    with Dataset(nc_path, 'r') as ds:
        for dim_name, dim in ds.dimensions.items():
            dim_sizes[dim_name] = len(dim)
        for vname, var in ds.variables.items():
            dims = tuple(var.dimensions)
            dtype = np.dtype(var.dtype.str)
            fill_value = getattr(var, '_FillValue', None)
            # Copy limited attrs (strings only) to avoid issues
            attrs: Dict[str, str] = {}
            for ak in var.ncattrs():
                try:
                    val = getattr(var, ak)
                    if isinstance(val, (str, int, float)):
                        attrs[ak] = str(val)
                except Exception:
                    pass
            var_meta[vname] = (dims, dtype, fill_value, attrs)
    return var_meta, dim_sizes


def select_variables(var_meta, selected: Optional[Iterable[str]]) -> Tuple[List[str], List[str]]:
    exclude_like = {'time', 'time_bnds', 'mcsec', 'date', 'datesec'}
    time_vars: List[str] = []
    static_vars: List[str] = []
    for name, (dims, dtype, _fill, _attrs) in var_meta.items():
        if name in exclude_like:
            continue
        if 'time' in dims and np.issubdtype(dtype, np.number):
            time_vars.append(name)
        elif 'time' not in dims:
            static_vars.append(name)
    if selected is not None:
        sel = set(selected)
        time_vars = [v for v in time_vars if v in sel]
        missing = sorted(sel.difference(set(time_vars).union(static_vars)))
        if missing:
            print(f"[warn] Selected variables not found or not time-varying: {missing}", file=sys.stderr)
    return time_vars, static_vars


def read_month_value(nc_path: str, var_name: str, dims: Tuple[str, ...], dtype: np.dtype, fill_value: Optional[float]) -> np.ndarray:
    with Dataset(nc_path, 'r') as ds:
        var = ds.variables[var_name]
        data = var[...]
        # Squeeze time dimension if present
        if 'time' in dims:
            time_axis = dims.index('time')
            data = np.take(data, indices=0, axis=time_axis)
        # Robust handling for NumPy 2.0 and masked arrays
        arr = np.asanyarray(data, dtype=np.float64)
        if np.ma.isMaskedArray(arr):
            arr = np.ma.filled(arr, np.nan)
        if fill_value is not None:
            try:
                fv = float(fill_value)
            except Exception:
                fv = fill_value
            arr = np.where(arr == fv, np.nan, arr)
        return arr


def ensure_output_dims(out_ds: Dataset, dim_sizes: Dict[str, int]) -> None:
    # Always create time dim length 1 for per-year files
    if 'time' not in out_ds.dimensions:
        out_ds.createDimension('time', 1)
    for dim, size in dim_sizes.items():
        if dim == 'time':
            continue
        if dim not in out_ds.dimensions:
            out_ds.createDimension(dim, size)


def write_static_vars(out_ds: Dataset, ref_path: str, static_vars: List[str], var_meta, compress_level: int) -> None:
    with Dataset(ref_path, 'r') as ref:
        for name in static_vars:
            dims, dtype, fill_value, attrs = var_meta[name]
            if name in out_ds.variables:
                continue
            zargs = {}
            v = out_ds.createVariable(name, dtype, dims, **zargs)
            # Copy minimal attrs
            for k, vval in attrs.items():
                if k in {'_FillValue'}:
                    continue
                try:
                    setattr(v, k, vval)
                except Exception:
                    pass
            # Write data
            v[...] = ref.variables[name][...]


def create_out_var(out_ds: Dataset, name: str, dims: Tuple[str, ...], src_dtype: np.dtype, fill_value: Optional[float], compress_level: int):
    # Insert 'time' dim if not already present (should be present in dims for time-varying vars)
    out_dims: Tuple[str, ...] = dims
    zargs = {}
    if fill_value is not None:
        zargs['fill_value'] = fill_value
    # Compression setup
    if compress_level and compress_level > 0:
        zargs.update({'zlib': True, 'complevel': int(compress_level), 'shuffle': True})
    v = out_ds.createVariable(name, src_dtype, out_dims, **zargs)
    return v


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Annual mean (no-leap weighted) using netCDF4 only, per-year outputs.")
    p.add_argument('inputs', nargs='+', help="Input glob patterns or file paths (e.g., '/path/case.elm.h0.????-??.nc').")
    p.add_argument('--vars', '--var', dest='variables', nargs='+', default=None, help='Variables to include (time-varying).')
    p.add_argument('--output', '-o', required=True, help='Output directory for per-year files.')
    p.add_argument('--compress-level', '-c', type=int, default=5, help='Compression level (0 disables). Default: 5.')
    p.add_argument('--allow-missing-months', action='store_true', help='Allow years with <12 months; weight by available months.')
    p.add_argument('--start-year', type=int, default=None, help='Inclusive start year to process (optional).')
    p.add_argument('--end-year', type=int, default=None, help='Inclusive end year to process (optional).')
    args = p.parse_args(argv)

    input_files = expand_input_patterns(args.inputs)
    if not input_files:
        print('[error] No input files matched.', file=sys.stderr)
        return 2

    year_map = group_files_by_year(input_files)
    if not year_map:
        print('[error] Could not group input files by year.', file=sys.stderr)
        return 3

    os.makedirs(args.output, exist_ok=True)

    # Validate year range if provided
    if args.start_year is not None and args.end_year is not None and args.start_year > args.end_year:
        print('[error] --start-year cannot be greater than --end-year.', file=sys.stderr)
        return 5

    # Use first file as reference for dimensions/vars
    ref_file = next(iter(sorted(input_files)))
    var_meta, dim_sizes = get_reference_metadata(ref_file)
    time_vars, static_vars = select_variables(var_meta, args.variables)
    if not time_vars:
        print('[error] No eligible time-varying numeric variables found.', file=sys.stderr)
        return 4

    case_stub = infer_case_stub(ref_file)

    # Filter years to requested range
    all_years = sorted(year_map.keys())
    years_to_process = [
        y for y in all_years
        if (args.start_year is None or y >= args.start_year) and (args.end_year is None or y <= args.end_year)
    ]
    if not years_to_process:
        print('[warn] No years to process after applying range filter.', file=sys.stderr)
        return 0

    for year in years_to_process:
        files = year_map[year]
        # Parse month index for weights
        months: List[int] = []
        for f in files:
            m = re.search(r"(\d{4})-(\d{2})\.nc$", os.path.basename(f))
            if m:
                months.append(int(m.group(2)))
            else:
                months.append(None)  # type: ignore
        if any(m is None for m in months):
            print(f"[warn] Skipping year {year} due to unparseable month in filenames.", file=sys.stderr)
            continue
        # Map month->file
        month_to_file = {m: f for m, f in zip(months, files)}
        present_months = sorted(month_to_file.keys())
        if len(present_months) < 12 and not args.allow_missing_months:
            print(f"[warn] Year {year} has {len(present_months)} months; skipping (use --allow-missing-months).", file=sys.stderr)
            continue

        # Prepare accumulators per variable
        accumulators: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        # Initialize shapes from first available file for each var
        sample_file = files[0]
        for vname in time_vars:
            dims, dtype, fill_value, _attrs = var_meta[vname]
            sample = read_month_value(sample_file, vname, dims, dtype, fill_value)
            num = np.zeros_like(sample, dtype=np.float64)
            den = np.zeros_like(sample, dtype=np.float64)
            accumulators[vname] = (num, den)

        # Accumulate weighted sums
        for month in present_months:
            f = month_to_file[month]
            weight = NOLEAP_DAYS_IN_MONTH[month - 1]
            for vname in time_vars:
                dims, dtype, fill_value, _attrs = var_meta[vname]
                val = read_month_value(f, vname, dims, dtype, fill_value)
                is_valid = np.isfinite(val)
                num, den = accumulators[vname]
                np.add(num, np.where(is_valid, val * weight, 0.0), out=num)
                np.add(den, np.where(is_valid, weight, 0.0), out=den)

        # Compute annual means
        annual_means: Dict[str, np.ndarray] = {}
        for vname, (num, den) in accumulators.items():
            with np.errstate(invalid='ignore', divide='ignore'):
                mean = np.divide(num, den, out=np.full_like(num, np.nan, dtype=np.float64), where=den > 0)
            # Cast back to source precision if float32, else keep float64
            _dims, src_dtype, fill_value, _attrs = var_meta[vname]
            if src_dtype == np.float32:
                mean = mean.astype(np.float32)
            annual_means[vname] = mean

        # Write per-year file
        if case_stub:
            out_name = f"{case_stub}.elm.h0.{year:04d}.annual_mean.nc"
        else:
            out_name = f"annual_mean_{year:04d}.nc"
        out_path = os.path.join(args.output, out_name)

        with Dataset(out_path, 'w', format='NETCDF4') as out_ds:
            ensure_output_dims(out_ds, dim_sizes)

            # time coordinate as integer year
            tvar = out_ds.createVariable('time', np.int32, ('time',))
            tvar[:] = np.array([year], dtype=np.int32)
            tvar.long_name = 'year'

            # Write static variables (coords/labels)
            write_static_vars(out_ds, ref_file, static_vars, var_meta, args.compress_level)

            # Write data variables
            for vname in time_vars:
                dims, src_dtype, fill_value, attrs = var_meta[vname]
                # Ensure 'time' is a dimension (it should be)
                if 'time' not in dims:
                    continue
                vout = create_out_var(out_ds, vname, dims, src_dtype, fill_value, args.compress_level)
                # Copy selected attrs
                for k, vval in attrs.items():
                    if k in {'_FillValue'}:
                        continue
                    try:
                        setattr(vout, k, vval)
                    except Exception:
                        pass
                # Write data with a time axis of length 1 at the correct position
                data = annual_means[vname]
                time_axis = dims.index('time')
                expanded = np.expand_dims(data, axis=time_axis)
                vout[...] = expanded

        print(f"[write] {out_path}")

    return 0


if __name__ == '__main__':
    raise SystemExit(main())


