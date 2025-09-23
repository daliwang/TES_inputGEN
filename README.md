### TESSFA — TES Input Generation Toolkit for ELM/E3SM

TESSFA provides a small, pragmatic toolkit to generate input datasets for E3SM/ELM single-site or regional simulations. It focuses on Area-Of-Interest (AOI) workflows for:

- Domain creation
- Surface datasets (surfdata)
- Meteorological forcing (ERA5, GSWP3)
- Utility scripts for linking, visualization, and QA/QC


## Repository layout

```
TES_inputGEN/
  TES_domainGEN.py                  # Domain generation entry point
  TES_surfdataGEN.py               # Legacy surfdata generator
  TES_surfdataGENv2.py             # Preferred surfdata generator (v2)
  TES_forcingGEN.py                # Generic forcing generator driver
  TES_ERA5forcingGEN.sub           # Slurm submission for ERA5 forcing workflow
  TES_GSWP3forcingGEN.sub          # Slurm submission for GSWP3 forcing workflow
  TES_AOI_forcingGEN.py            # AOI-specific forcing pipeline
  TES_AOI_testbench.ipynb          # AOI testing/verification notebook
  TES_surfdataGENv2.ipynb          # Notebook for surfdata v2 workflow
  ncdiff2.py                       # NetCDF helper/utility
  TES_AOI_scripts/                 # AOI-focused helper scripts
    domain_creation.sh
    forcing_creation.sbatch
    forcing_domain_link_creation.py
    forcinglink_creation.py
    shape2gridID.py
    Show2DLocation.py
    Show2DVariable.py
    Show2DVariables.v2.py
    Softlink_creation.py
    TES_AOI_domainGEN.py
    TES_AOI_forcingGEN.py
    TES_AOI_surfdataGEN.py
    TES_TNgridID.py
    Varible2Geotiff.py
```


## Quickstart

1) Install system tools (macOS/Homebrew shown)

```bash
brew update
brew install nco cdo gdal proj
```

2) Create a conda environment (recommended)

```bash
conda create -n tessfa -c conda-forge \
  python=3.10 numpy pandas xarray netcdf4 scipy dask rasterio rioxarray \
  geopandas shapely pyproj fiona gdal matplotlib cartopy
conda activate tessfa
```

3) Verify scripts are available and see usage

```bash
cd TES_inputGEN
python TES_domainGEN.py --help
python TES_surfdataGENv2.py --help
python TES_forcingGEN.py --help
```

4) Recommended inputs you will need

- AOI polygon (ESRI Shapefile or GeoJSON) in a projected CRS
- Access to land surface data sources used by surfdata generation
- Access to ERA5 or GSWP3 meteorological datasets (local paths or remote)


## Typical workflows

### 1) Create a domain for your AOI

- Use `TES_inputGEN/TES_domainGEN.py` for general domain generation, or the AOI-specific driver `TES_inputGEN/TES_AOI_scripts/TES_AOI_domainGEN.py`.
- These scripts typically take an AOI vector file, a target grid definition, and output an ELM/E3SM-compatible domain NetCDF.

Example (pattern — use `--help` to see actual flags):

```bash
python TES_inputGEN/TES_domainGEN.py \
  --aoi /path/to/aoi.shp \
  --dx 0.125 --dy 0.125 \
  --out /path/to/output/domain.nc
```


### 2) Generate surfdata for the AOI

- Prefer `TES_inputGEN/TES_surfdataGENv2.py` (v2), with the legacy `TES_surfdataGEN.py` as a fallback.
- Provide the domain file and any catalog/configs the script expects.

Example (pattern — consult `--help`):

```bash
python TES_inputGEN/TES_surfdataGENv2.py \
  --domain /path/to/output/domain.nc \
  --out /path/to/output/surfdata.nc
```


### 3) Generate meteorological forcing

There are two common sources supported in this repo: ERA5 and GSWP3.

- Generic driver: `TES_inputGEN/TES_forcingGEN.py`
- AOI driver: `TES_inputGEN/TES_AOI_forcingGEN.py`
- Slurm submissions: `TES_inputGEN/TES_ERA5forcingGEN.sub`, `TES_inputGEN/TES_GSWP3forcingGEN.sub`

Run locally (pattern):

```bash
python TES_inputGEN/TES_forcingGEN.py \
  --domain /path/to/output/domain.nc \
  --source ERA5 \
  --start-year 2000 --end-year 2010 \
  --outdir /path/to/output/forcing
```

Run via Slurm (pattern):

```bash
sbatch TES_inputGEN/TES_ERA5forcingGEN.sub
sbatch TES_inputGEN/TES_GSWP3forcingGEN.sub
```

Use `--help` on the Python drivers and open the `.sub` files to set account/queue, years, paths, and modules for your cluster.


### 4) Link, visualize, and QA/QC

- Link forcing and domain: `TES_inputGEN/TES_AOI_scripts/forcing_domain_link_creation.py`, `forcinglink_creation.py`, `Softlink_creation.py`
- Map AOI to grid IDs/tiles: `shape2gridID.py`, `TES_TNgridID.py`
- Visualization: `Show2DLocation.py`, `Show2DVariable.py`, `Show2DVariables.v2.py`
- Export to GeoTIFF: `Varible2Geotiff.py`

Most of these scripts support `--help`. They are useful for verifying alignment, checking variable values, and producing quicklooks.


## Notebooks

- `TES_inputGEN/TES_AOI_testbench.ipynb`: End-to-end AOI sanity checks and diagnostics
- `TES_inputGEN/TES_surfdataGENv2.ipynb`: Interactive surfdata v2 workflow exploration

Launch with Jupyter:

```bash
conda activate tessfa
jupyter lab
```


## Installation notes (macOS)

- Install `gdal` before Python packages that depend on it. The conda-forge stack above ensures compatible versions of `gdal`, `geopandas`, `fiona`, and `rasterio`.
- On Apple Silicon, prefer conda-forge; avoid mixing with system Python or pip wheels compiled against a different `gdal`/`proj`.
- For NetCDF command-line work, `nco` and `cdo` are installed via Homebrew above.


## Troubleshooting

- If you see `PROJ`/`GDAL` errors, ensure the conda environment is active and that `gdal`/`proj` come from conda-forge.
- If NetCDF writes fail, verify `netcdf4` is installed and that output directories exist and are writable.
- For Slurm runs, edit `.sub` files to match your cluster's modules, account, QoS/partition, and filesystem paths.


## Script reference (brief)

- `TES_domainGEN.py`: Create ELM/E3SM domain NetCDF from AOI and grid settings.
- `TES_surfdataGENv2.py`: Generate surfdata; use this version when possible.
- `TES_surfdataGEN.py`: Legacy surfdata generator retained for compatibility.
- `TES_forcingGEN.py`: Driver for building forcing from supported sources.
- `TES_AOI_forcingGEN.py`: AOI-oriented forcing workflow glue code.
- `TES_ERA5forcingGEN.sub`, `TES_GSWP3forcingGEN.sub`: Slurm submit scripts for ERA5/GSWP3.
- `ncdiff2.py`: NetCDF utility script (e.g., variable diffs/merges/helpers).
- `TES_AOI_scripts/`: Helpers for domain creation, linking, visualization, and exports.


## Contributing

- Please add `--help`/argparse usage to new scripts and keep AOI examples up to date.
- Prefer `conda-forge` packages and pin major versions in notebooks for reproducibility.
- Open PRs with a short description of the workflow addressed and a minimal example.


## License

Specify your project's license here (e.g., MIT, BSD-3-Clause). If unsure, add a `LICENSE` file in the repo root and update this section.


