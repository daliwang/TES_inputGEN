## Annual Mean Computation (ELM outputs)

This folder contains two scripts to compute annual means from ELM NetCDF outputs.
- **compute_annual_means_nc4.py**: netCDF4-only, streams monthly files per year (low memory), writes one file per year.
- **compute_annual_means.py**: xarray-based, supports single-file or multi-file inputs, optional dask chunking.

### Requirements
- For `compute_annual_means_nc4.py`: numpy, netCDF4
- For `compute_annual_means.py`: numpy, xarray (optionally dask, and a backend like netCDF4 or h5netcdf)

Tip: When using zsh, quote globs so the script receives the pattern.

### Quick start — write all annual means next to the monthly files
Use the netCDF4 streaming script. It groups monthly h0 files by year and writes `<case>.elm.h0.YYYY.annual_mean.nc` in the same directory.

```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/compute_annual_means_nc4.py \
  "/path/to/run_dir/*.elm.h0.????-??.nc" \
  --output "/path/to/run_dir" -c 5
```

- Replace `/path/to/run_dir` with the directory containing monthly `*.elm.h0.YYYY-MM.nc` files.
- By default it averages all time-varying numeric variables.

Optional flags:
- `--vars NAME1 NAME2 ...`: limit to specific variables (space-separated)
- `--allow-missing-months`: process years with <12 months (no-leap weighted)
- `-c/--compress-level N`: compression level (0 disables; default 5)

Example with options:
```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/compute_annual_means_nc4.py \
  "/path/to/run_dir/*.elm.h0.????-??.nc" \
  --vars TBOT FSDS GPP FLDS RH2M \
  --output "/path/to/run_dir" \
  -c 5 --allow-missing-months
```

### xarray workflow — single file or specific monthly set
`compute_annual_means.py` computes a mean over the `time` dimension for selected variables and writes a single output file.

- Single-file mode (h1 or h0 time series):
```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/compute_annual_means.py \
  --in /path/to/input_time_series.nc \
  --out /path/to/output_annual.nc \
  --vars GPP,NPP,TBOT \
  --chunks auto
```

- Monthly multi-file mode for one year (glob or explicit list):
```bash
# Using a glob for one year's h0 files
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/compute_annual_means.py \
  --monthly-glob "/path/to/run_dir/case.elm.h0.0005-*.nc" \
  --out "/path/to/run_dir/case.elm.h0.0005_annual.nc" \
  --vars GPP,NPP,TBOT

# Or list files explicitly (comma-separated)
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/compute_annual_means.py \
  --monthly-files "/path/run/case.elm.h0.0005-01.nc,/path/run/case.elm.h0.0005-02.nc,..." \
  --out "/path/to/run_dir/case.elm.h0.0005_annual.nc"
```

Common options:
- `--vars name1,name2,...`: comma-separated variables to include (default: all numeric with `time`)
- `--drop-time-dim`: drop the `time` dimension in output (single snapshot)
- `--chunks auto` (default) or `--chunks time:168,lat:100,lon:100`; use `--no-dask` to disable
- `--engine netcdf4` or `--engine h5netcdf` to select backend

### Output naming
- `compute_annual_means_nc4.py`: writes one file per year in the output directory (same dir recommended), named like `<case>.elm.h0.YYYY.annual_mean.nc`.
- `compute_annual_means.py`: writes a single file; if `--out` is omitted, it derives a name based on the input (e.g., `..._annual.nc`).

### Notes
- nc4 script weights months by no-leap day counts; xarray script computes simple mean over `time` and adds `cell_methods: time: mean (annual)`.
- Always quote globs: "`/path/*.elm.h0.????-??.nc`".

### Batch plotting annual-mean variables
Use `batch_annual_plots.py` to generate plots for all time-varying variables from each `*.annual_mean.nc` file in a directory. It auto-detects variables that include a `time` dimension and calls `Show2DVariablesBatch.py`.

```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/batch_annual_plots.py \
  --run-dir /path/to/run_dir \
  --mask /path/to/mask.nc \
  --out /path/to/plots \
  --time-index 0 \
  --cmap viridis \
  --dpi 150 \
  --year-digits 4
```

- `--run-dir`: directory containing `*.annual_mean.nc` files
- `--mask`: mask NetCDF with variables `x`, `y`, `gridID`
- `--out`: output directory for plots (created if missing)
- Optional: `--cmap`, `--dpi`, `--year-digits`, and `--time-index` (annual files typically have length-1 time)

The script invokes `Show2DVariablesBatch.py` at:
`/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_AOI_scripts/Show2DVariablesBatch.py`.
