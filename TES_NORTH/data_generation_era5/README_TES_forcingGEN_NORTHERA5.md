## TES_forcingGEN_NORTHERA5

Generate 1D TES forcing from Daymet+ERA5 TESSFA North inputs. Processes NetCDF forcing in nested subdirectories and writes 1D outputs with standardized names.

### What it does
- Reads 2D forcing files and extracts values over land cells into a compact 1D layout (time, ni), plus metadata variables (gridID, LATIXY, LONGXY, time).
- Walks `input_path` recursively, ignoring files directly in the top-level of `input_path`. Only files ending with `.nc` are processed.
- Mirrors the input folder structure under `output_path` (per-variable and other subfolders preserved).
- Removes leap day if a 232-step sequence is detected (232 -> 224) for consistency.
- Builds output filenames in the form: `climforc.<DATASETID>.<resolution>.1d.<variable>.<yyyy-mm>.nc`.

### File naming
- Input filename pattern expected (source):
  - `clmforc.<dataset>.<resolution_km>.<var>.<yyyy-mm>.nc`
  - Example: `clmforc.Daymet_ERA5_TESSFA.4km.Prec.2023-04.nc`
- Output filename pattern (target):
  - `climforc.<DATASETID>.<resolution>.1d.<var>.<yyyy-mm>.nc`
  - `<DATASETID>` comes from the optional CLI arg or is derived from `output_path` when omitted.
  - `<resolution>` is inferred from the input filename token ending with `km` (fallback `4km`).

### Dependencies
- Python 3.11+
- Packages: `netCDF4`, `numpy`

### Python script usage
File: `TES_forcingGEN_NORTHERA5.py`

```
python3 TES_forcingGEN_NORTHERA5.py <input_path> <output_path> <time_steps> [DATASETID]
```
- `<input_path>`: root folder containing nested subfolders with `.nc` inputs. Top-level files are ignored.
- `<output_path>`: destination root; subfolder structure mirrors `input_path`.
- `<time_steps>`: integer number of timesteps to write per file (`-1` for full time series).
- `[DATASETID]` (optional): dataset id used in output filenames (e.g., `Daymet_ERA5_TESSFA_NORTH`). If omitted, it is derived from `output_path` (expects `.../<DATASETID>/entire_domain/forcing/`).

Examples:
```
# Process a single timestep, local paths
python3 TES_forcingGEN_NORTHERA5.py ./temp ./forcing 1 Daymet_ERA5_TESSFA_NORTH

# Process all timesteps (no limit)
python3 TES_forcingGEN_NORTHERA5.py /path/to/input /path/to/output -1 Daymet_ERA5_TESSFA_NORTH
```

Notes:
- Outputs are written under `output_path/<relative_subdir_from_input>/climforc.<DATASETID>.<resolution>.1d.<var>.<yyyy-mm>.nc`.
- If the input file has 232 time steps, the script normalizes to 224 (leap-day removal). The `time` coordinate is written with units `days since <year>-<month>-01 00:00:00` and values correspond to day-of-month with 8 sub-daily intervals.

### Slurm wrapper usage
File: `TES_forcingGEN_NORTHERA5.sub`

- Loads environment, sets variables, cleans `forcing_output_path`, and invokes the Python generator.
- Relevant variables to edit:
  - `time_steps` (e.g., `-1` for all, or `1` for a quick test)
  - `DATASETID` (e.g., `Daymet_ERA5_TESSFA_NORTH`)
  - `TES_forcing_path` (input root)
  - `forcing_output_path` (output root; can be absolute or `./forcing/` for local testing)

Local run (no Slurm):
```
sh TES_forcingGEN_NORTHERA5.sub
```

Slurm submission (if using srun in the script):
```
sbatch TES_forcingGEN_NORTHERA5.sub
```

### Troubleshooting
- Only one output file appears in `output_path`:
  - Check subdirectories. Outputs mirror the input directory structure and are written inside corresponding subfolders.
- Filename doesn’t include the correct dataset id:
  - Pass `DATASETID` as the 4th CLI argument, or ensure `output_path` matches pattern `.../<DATASETID>/entire_domain/forcing/` for automatic derivation.
- IndexError: size of data array does not conform to slice:
  - Ensure `<time_steps>` matches your intent. The script slices the time vector to `time_steps` before writing.
