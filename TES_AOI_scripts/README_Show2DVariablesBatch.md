## Show2DVariablesBatch.py — Plot 2D variables (PNG) and export GeoTIFFs

Batch plot variables from a NetCDF "data" file onto a 2D grid defined by a NetCDF "mask" (with `x`, `y`, `gridID`). The script maps 1D active-cell values back to a 2D grid using `gridID`, supports PNG exports, optional variable×area visualization, and multi-band GeoTIFFs.

### Requirements
- Python 3, `numpy`, `netCDF4`, `matplotlib` (non-interactive backend)
- Optional for GeoTIFF export: `rasterio`

### Inputs
- Mask NetCDF must contain: `x`, `y`, `gridID` (gridID indexes flattened 2D grid)
- Data NetCDF must contain the variables to plot (1D over active cells or flattenable to grid), optional `time` dimension
- For GeoTIFF export: data should provide `lon` and `lat` (or equivalents) to georeference

### Usage
```bash
python3 Show2DVariablesBatch.py \
  --mask /path/to/mask.nc \
  --data /path/to/data.nc \
  --vars VAR1 VAR2 ... \
  --output-dir /path/to/plots \
  [--time-index 0] [--prefix PREFIX] [--cmap viridis] [--dpi 150] \
  [--vmin MIN] [--vmax MAX] \
  [--x-min 1 --x-max XMAX --y-min 1 --y-max YMAX] \
  [--year-digits 4] [--year YYYY]
```

### PNG: variable × area for visualization
Multiply plotted values by gridcell area for PNGs (e.g., `kgC m-2 × m2` = `kgC per cell`).
```bash
python3 Show2DVariablesBatch.py \
  --mask /path/to/mask.nc \
  --data /path/to/carbon_loss_0800.nc \
  --vars carbon_loss \
  --output-dir /path/to/plots \
  --png-multiply-by-area \
  --time-index 0
```
- Colorbar label shows units (e.g., `carbon_loss x area [kgC m-2 * m2]`).
- Output filename includes `_xarea` when applied.

### GeoTIFF export (multi-band)
Export raw variable and optional extra bands for `area` and `variable×area`.
```bash
python3 Show2DVariablesBatch.py \
  --mask /path/to/mask.nc \
  --data /path/to/carbon_loss_0800.nc \
  --vars carbon_loss \
  --export-geotiff --tiff-dir /path/to/geotiff \
  --tiff-include-area --tiff-include-product \
  --time-index 0
```
GeoTIFF bands:
- band 1: raw variable (e.g., `carbon_loss`)
- band 2: `area` (if `--tiff-include-area`)
- band 3: `variable × area` (if `--tiff-include-product`)

### Subsetting and appearance
- `--x-min/--x-max/--y-min/--y-max`: 1-based inclusive indices for subdomain plotting
- `--cmap`, `--vmin`, `--vmax`, `--dpi`: control appearance

### Year tagging in filenames
Filenames are tagged using, in order of priority:
1) `--year` if provided
2) `time` coordinate from data (at `--time-index`)
3) last 4–5 digit sequence detected in the data filename

### Notes
- Units: colorbar labels include units retrieved from the data. When using `--png-multiply-by-area`, the label composes both units.
- GeoTIFF CRS: exports as `EPSG:4326` based on `lon`/`lat` from the data file.
- Quote paths/globs in zsh.

