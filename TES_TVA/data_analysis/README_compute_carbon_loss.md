## Compute and Plot Annual Carbon Loss

This documents `compute_carbon_loss.py` and how to visualize/export the results using `Show2DVariablesBatch.py`.

### 1) Compute carbon_loss from annual_mean NetCDF
Formula: `carbon_loss = WOODC + CWDC + TOTLITC + 0.1 * TOTSOILC`.

```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_analysis/compute_carbon_loss.py \
  --annual /path/to/run/annual/case.elm.h0.0800.annual_mean.nc \
  --mask   /path/to/mask.nc \
  --out    /path/to/output/carbon_loss_0800.nc
```

Output (1D over active cells, aligned to `gridID`):
- `carbon_loss` (units propagated from source pools; often `kgC m-2`)
- `gridID` (index of active cells)
- `lon`, `lat` (degrees)
- `area` (m2)

### 2) Plot PNGs (variable × area for visualization)
Use `--png-multiply-by-area` to show totals per gridcell (e.g., kgC per cell):

```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_AOI_scripts/Show2DVariablesBatch.py \
  --mask /path/to/mask.nc \
  --data /path/to/output/carbon_loss_0800.nc \
  --vars carbon_loss \
  --output-dir /path/to/plots \
  --png-multiply-by-area \
  --time-index 0
```

### 3) Export GeoTIFFs (multi-band)
Export raw `carbon_loss` plus `area`, and optionally the product as an additional band:

```bash
python3 /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_AOI_scripts/Show2DVariablesBatch.py \
  --mask /path/to/mask.nc \
  --data /path/to/output/carbon_loss_0800.nc \
  --vars carbon_loss \
  --export-geotiff --tiff-dir /path/to/geotiff \
  --tiff-include-area --tiff-include-product \
  --time-index 0
```

GeoTIFF bands:
- band 1: `carbon_loss` (raw)
- band 2: `area` (if `--tiff-include-area`)
- band 3: `carbon_loss × area` (if `--tiff-include-product`)

### Notes
- The masking/grid mapping assumes `gridID` indexes the flattened 2D grid.
- Quote paths/globs in zsh.

