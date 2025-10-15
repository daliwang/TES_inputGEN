#!/usr/bin/env python3
import argparse
import datetime as dt
import os
from typing import Optional, Tuple, List

import numpy as np
from netCDF4 import Dataset


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Compute annual carbon_loss from annual_mean NetCDF and write 1D arrays aligned to mask gridIDs.')
    p.add_argument('--annual', required=True, help='Path to annual_mean NetCDF (e.g., case.elm.h0.0800.annual_mean.nc)')
    p.add_argument('--mask', required=True, help='Path to mask NetCDF (must contain x, y, gridID) to define grid shape and active indices')
    p.add_argument('--out', required=True, help='Output NetCDF to write carbon_loss, area, gridID, lon, lat (1D over active cells)')
    p.add_argument('--time-index', type=int, default=0, help='Time index in annual file (usually 0). Default: 0')
    return p.parse_args()


def load_mask(mask_path: str) -> Tuple[np.ndarray, Tuple[int, int]]:
    ds = Dataset(mask_path)
    try:
        x = ds.variables['x'][:]
        y = ds.variables['y'][:]
        gridID = ds.variables['gridID'][:]
        if hasattr(gridID, 'compressed'):
            gridID = gridID.compressed()
        gridID = np.asarray(gridID, dtype=np.int64)
        shape2d = (len(y), len(x))
        return gridID, shape2d
    finally:
        ds.close()


def get_units(ds: Dataset, var_name: str) -> Optional[str]:
    if var_name not in ds.variables:
        return None
    v = ds.variables[var_name]
    for key in ('units', 'Units', 'unit', 'UNIT'):
        try:
            u = getattr(v, key)
            if isinstance(u, str) and u.strip():
                return u.strip()
        except Exception:
            pass
    return None


def read_var_aligned_1d(ds: Dataset, var_name: str, time_index: int, gridID: np.ndarray, shape2d: Tuple[int, int]) -> Optional[np.ndarray]:
    if var_name not in ds.variables:
        return None
    v = ds.variables[var_name]
    data = v[:]
    if hasattr(v, 'dimensions') and 'time' in v.dimensions:
        t_axis = v.dimensions.index('time')
        data = np.take(data, indices=time_index, axis=t_axis)
    arr = np.asarray(data)
    arr = np.squeeze(arr)
    # Case A: already 1D per active cell
    if arr.ndim == 1 and arr.size == gridID.size:
        return arr.astype(np.float64)
    # Case B: 1D flattened full grid (size == Ny*Nx) -> select by gridID
    if arr.ndim == 1 and arr.size == (shape2d[0] * shape2d[1]):
        return arr.astype(np.float64)[gridID]
    # Case C: 2D grid -> flatten then select by gridID
    if arr.ndim == 2 and arr.shape == shape2d:
        return np.ravel(arr.astype(np.float64))[gridID]
    return None


def pick_first_present(ds: Dataset, names: List[str], time_index: int) -> Optional[np.ndarray]:
    for n in names:
        v = read_var_aligned_1d(ds, n, time_index, np.arange(0), (0, 0))
        if v is not None:
            return v
    return None


def main() -> int:
    args = parse_args()

    gridID, shape2d = load_mask(args.mask)

    with Dataset(args.annual) as ds:
        woodc = read_var_aligned_1d(ds, 'WOODC', args.time_index, gridID, shape2d)
        cwdc = read_var_aligned_1d(ds, 'CWDC', args.time_index, gridID, shape2d)
        totlitc = read_var_aligned_1d(ds, 'TOTLITC', args.time_index, gridID, shape2d)
        totcolc = read_var_aligned_1d(ds, 'TOTCOLC', args.time_index, gridID, shape2d)

        if any(a is None for a in (woodc, cwdc, totlitc, totcolc)):
            missing = [n for n, a in [('WOODC', woodc), ('CWDC', cwdc), ('TOTLITC', totlitc), ('TOTCOLC', totcolc)] if a is None]
            raise ValueError(f"Missing required variables in annual file or unexpected shapes: {', '.join(missing)}")

        carbon_loss_1d = woodc + cwdc + totlitc + (0.1 * totcolc)

        # Coordinates/area as 1D aligned to gridID
        lon_1d = None
        lat_1d = None
        area_1d = None
        def first_var_aligned_1d(name_list: List[str]) -> Optional[np.ndarray]:
            for n in name_list:
                arr = read_var_aligned_1d(ds, n, args.time_index, gridID, shape2d)
                if arr is not None:
                    return arr
            return None

        lon_1d = first_var_aligned_1d(['lon', 'LONGXY', 'LON'])
        lat_1d = first_var_aligned_1d(['lat', 'LATIXY', 'LAT'])
        area_1d = first_var_aligned_1d(['area', 'AREA', 'AREA_GRID'])

        if lon_1d is None:
            lon_1d = np.full_like(carbon_loss_1d, np.nan, dtype=np.float64)
        if lat_1d is None:
            lat_1d = np.full_like(carbon_loss_1d, np.nan, dtype=np.float64)
        if area_1d is None:
            area_1d = np.full_like(carbon_loss_1d, np.nan, dtype=np.float64)

        # Determine output units from source pools (assume same units across pools)
        pool_units = [
            get_units(ds, 'WOODC'),
            get_units(ds, 'CWDC'),
            get_units(ds, 'TOTLITC'),
            get_units(ds, 'TOTCOLC'),
        ]
        non_empty = [u for u in pool_units if u]
        if non_empty and len(set(non_empty)) == 1:
            out_units = non_empty[0]
        elif non_empty:
            # Units disagree; fall back to first non-empty
            out_units = non_empty[0]
        else:
            out_units = 'unknown'

    # Derive output filename to include annual file info when a directory or stem is provided
    annual_base = os.path.basename(args.annual)
    annual_stem = os.path.splitext(annual_base)[0]
    derived_name = f"{annual_stem}_carbon_loss.nc"

    if os.path.isdir(args.out) or str(args.out).endswith(os.sep):
        out_path = os.path.join(args.out, derived_name)
    elif os.path.splitext(str(args.out))[1].lower() == '.nc':
        out_path = str(args.out)
    else:
        out_dir = os.path.dirname(str(args.out))
        if out_dir:
            out_path = os.path.join(out_dir, derived_name)
        else:
            out_path = derived_name

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or '.', exist_ok=True)
    with Dataset(out_path, 'w', format='NETCDF4') as out:
        # Copy global attributes from source and annotate
        try:
            for ak in ds.ncattrs():
                try:
                    out.setncattr(ak, getattr(ds, ak))
                except Exception:
                    pass
        except Exception:
            pass

        title_prefix = 'Annual carbon loss (derived)'
        try:
            existing_title = getattr(out, 'title', None)
            out.title = f"{title_prefix} — {existing_title}" if existing_title else title_prefix
        except Exception:
            pass

        try:
            ts = dt.datetime.utcnow().isoformat(timespec='seconds') + 'Z'
            hist_line = f"{ts}: carbon_loss computed as WOODC + CWDC + TOTLITC + 0.1*TOTSOILC from {os.path.basename(args.annual)}"
            existing_history = getattr(out, 'history', None)
            out.history = (existing_history + '\n' + hist_line) if existing_history else hist_line
        except Exception:
            pass
        n = carbon_loss_1d.size
        out.createDimension('grid', n)

        v_grid = out.createVariable('gridID', np.int64, ('grid',))
        v_grid[:] = gridID

        v_cl = out.createVariable('carbon_loss', np.float64, ('grid',))
        v_cl.long_name = 'Annual carbon loss'
        v_cl.units = out_units
        v_cl[:] = carbon_loss_1d

        v_lon = out.createVariable('lon', np.float64, ('grid',))
        v_lon.units = 'degrees_east'
        v_lon[:] = lon_1d

        v_lat = out.createVariable('lat', np.float64, ('grid',))
        v_lat.units = 'degrees_north'
        v_lat[:] = lat_1d

        v_area = out.createVariable('area', np.float64, ('grid',))
        v_area.units = 'm2'
        v_area[:] = area_1d

    print(f"[write] {out_path}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())


