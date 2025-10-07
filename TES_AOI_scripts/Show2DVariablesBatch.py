#!/usr/bin/env python3
import argparse
import os
import re

import numpy as np
import netCDF4 as nc

# Use non-interactive backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description='Batch plot 2D variables to files using mask/data NetCDFs.')
    p.add_argument('--mask', required=True, help='Path to mask NetCDF (contains x, y, gridID).')
    p.add_argument('--data', required=True, help='Path to data NetCDF (variables to plot).')
    p.add_argument('--vars', nargs='+', required=True, help='List of variable names to plot.')
    p.add_argument('--time-index', type=int, default=0, help='Time index to plot when variable has time dimension. Default: 0')
    p.add_argument('--output-dir', '-o', required=True, help='Directory to write plot files.')
    p.add_argument('--prefix', default='', help='Optional filename prefix.')
    p.add_argument('--cmap', default='viridis', help='Matplotlib colormap. Default: viridis')
    p.add_argument('--dpi', type=int, default=150, help='Output DPI. Default: 150')
    p.add_argument('--vmin', type=float, default=None, help='Global vmin for color scale (optional).')
    p.add_argument('--vmax', type=float, default=None, help='Global vmax for color scale (optional).')
    p.add_argument('--x-min', type=int, default=1, help='1-based x start for subdomain. Default: 1')
    p.add_argument('--x-max', type=int, default=None, help='1-based x end for subdomain (inclusive). Default: full')
    p.add_argument('--y-min', type=int, default=1, help='1-based y start for subdomain. Default: 1')
    p.add_argument('--y-max', type=int, default=None, help='1-based y end for subdomain (inclusive). Default: full')
    p.add_argument('--year-digits', type=int, choices=[4, 5], default=4, help='Digits for year in filename (4 or 5). Default: 4')
    p.add_argument('--year', type=int, default=None, help='Explicit year override for filenames.')
    return p.parse_args()


def load_mask(mask_path: str):
    ds = nc.Dataset(mask_path)
    x_coords = ds.variables['x'][:]
    y_coords = ds.variables['y'][:]
    gridID = ds.variables['gridID'][:]
    # If masked, compress to 1D active indices
    if hasattr(gridID, 'compressed'):
        gridID = gridID.compressed()
    gridID = np.asarray(gridID, dtype=int)
    active_array = np.zeros((len(y_coords), len(x_coords)), dtype=int)
    flattened_active_array = active_array.flatten()
    # Mark active cells
    for gid in gridID:
        if 0 <= gid < flattened_active_array.size:
            flattened_active_array[gid] = 1
    reshaped_active_array = flattened_active_array.reshape(active_array.shape)
    return ds, x_coords, y_coords, gridID, reshaped_active_array


def read_variable(data_ds: nc.Dataset, var_name: str, time_index: int) -> np.ndarray:
    var = data_ds.variables[var_name]
    data = var[:]
    # Handle time dimension if present (assume time is first axis if present)
    if data.ndim == 3:
        if time_index < 0 or time_index >= data.shape[0]:
            raise IndexError(f'time_index {time_index} out of range for variable {var_name}')
        data = data[time_index]
    # If data has a singleton dimension (e.g., j=1), squeeze it
    data = np.squeeze(data)
    # If masked array, convert to 1D values of active grid cells
    if hasattr(data, 'compressed'):
        data1d = data.compressed()
    else:
        data1d = np.asarray(data)
    return data1d


def extract_year_from_filename(path: str) -> int | None:
    base = os.path.basename(path)
    # Prefer final .YYYY or .YYYYY before extension
    m = re.search(r"\.(\d{4,5})(?:\.|_)?[^.]*$", base)
    if not m:
        # Fallback to last 4-5 digit sequence anywhere in basename
        candidates = re.findall(r'(\d{4,5})', base)
        if candidates:
            mval = candidates[-1]
            try:
                return int(mval)
            except Exception:
                return None
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None
    return None


def map_to_2d(gridID: np.ndarray, values1d: np.ndarray, shape2d: tuple) -> np.ndarray:
    out = np.full(shape2d, np.nan, dtype=float)
    n = min(len(gridID), values1d.size)
    # Map first n active points
    for i in range(n):
        gid = gridID[i]
        if 0 <= gid < out.size:
            out[np.unravel_index(gid, shape2d)] = float(values1d[i])
    return out


def main() -> int:
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    mask_ds, x_coords, y_coords, gridID, active_mask_2d = load_mask(args.mask)
    try:
        data_ds = nc.Dataset(args.data)
        try:
            # Determine subdomain extents
            x_min = args.x_min
            x_max = args.x_max or len(x_coords)
            y_min = args.y_min
            y_max = args.y_max or len(y_coords)
            # Convert to 0-based slice indexes
            x0, x1 = x_min - 1, x_max
            y0, y1 = y_min - 1, y_max
            sub_x = x_coords[x0:x1]
            sub_y = y_coords[y0:y1]
            shape2d = active_mask_2d.shape

            # Determine year for naming (priority: --year > dataset time > filename > tNNN)
            year_val = args.year
            if year_val is None:
                try:
                    if 'time' in data_ds.variables:
                        t = data_ds.variables['time'][:]
                        if t.ndim >= 1 and t.shape[0] > args.time_index:
                            tv = t[args.time_index]
                            year_from_time = int(np.asarray(tv).item())
                            # Guard against zero/negative or obviously invalid
                            if year_from_time > 0:
                                year_val = year_from_time
                except Exception:
                    year_val = None
            if year_val is None:
                year_val = extract_year_from_filename(args.data)
            if year_val is not None and year_val > 0:
                year_full = f"{year_val:05d}"
                year_tag = year_full[-4:] if args.year_digits == 4 else year_full
            else:
                year_tag = f"t{args.time_index:03d}"

            for var_name in args.vars:
                if var_name not in data_ds.variables:
                    print(f"[warn] Variable not found: {var_name}")
                    continue
                try:
                    values1d = read_variable(data_ds, var_name, args.time_index)
                except Exception as e:
                    print(f"[warn] Skipping {var_name}: {e}")
                    continue

                # Map to 2D domain
                data2d = map_to_2d(gridID, values1d, shape2d)
                sub_data = data2d[y0:y1, x0:x1]

                plt.figure(figsize=(8, 6), dpi=args.dpi)
                mesh = plt.pcolormesh(sub_x, sub_y, sub_data, shading='auto', cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
                plt.colorbar(mesh, label=var_name)
                plt.title(f'{var_name} (t={args.time_index})')
                plt.xlabel('X Coordinates')
                plt.ylabel('Y Coordinates')

                prefix = (args.prefix + '_') if args.prefix else ''
                out_name = f"{prefix}{var_name}_{year_tag}.png"
                out_path = os.path.join(args.output_dir, out_name)
                plt.tight_layout()
                plt.savefig(out_path, dpi=args.dpi)
                plt.close()
                print(f"[write] {out_path}")
        finally:
            data_ds.close()
    finally:
        mask_ds.close()

    return 0


if __name__ == '__main__':
    raise SystemExit(main())


