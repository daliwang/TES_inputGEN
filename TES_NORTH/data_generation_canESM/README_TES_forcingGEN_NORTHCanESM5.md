## CanESM5 TES forcing generation

This folder contains per-case scripts to convert 2D CanESM5 forcing to 1D TES forcing for the NORTH domain, plus a Slurm batch script to run test and production workflows.

### Inputs

- Base input directory: `/gpfs/wolf2/cades/cli185/proj-shared/wangd/canESM5`
- Each case is a subfolder that starts with `CanESM5*` and contains:
  - `TPHWL3Hrly/`
  - `Solar3Hrly/`
  - `Precip3Hrly/`

### Scripts (per case)

- `TES_forcingGEN_NORTHCanESM5_1pctCO2-bgc_r1i1p1f1_DBCCA_Daymet_TESSFA1.py`
- `TES_forcingGEN_NORTHCanESM5_historical_r1i1p1f1_DBCCA_Daymet_TESSFA1.py`
- `TES_forcingGEN_NORTHCanESM5_piControl_r1i1p1f1_DBCCA_Daymet_TESSFA1.py`
- `TES_forcingGEN_NORTHCanESM5_ssp585_r1i1p1f1_DBCCA_Daymet_TESSFA1.py`

Each script writes files as:
`clmforc.<DATASETID>.<resolution>.1d.<var>.<period>.nc`

Optional year filter:
```
python <script> <input_path> <output_path> <time_steps> [start_year end_year]
```

### Batch script

`TES_NORTHCanESM5forcingGEN.sub` supports test and production modes.

**Quick test (1 case, 20 years, 1 timestep per file):**
```
sbatch TES_NORTHCanESM5forcingGEN.sub
```

Optional overrides for test:
```
sbatch --export=MODE=test,TEST_CASE_INDEX=2,TEST_START_YEAR=1980,TEST_YEAR_COUNT=20 \
  TES_NORTHCanESM5forcingGEN.sub
```

**Production (4 cases in parallel, all timesteps):**
```
sbatch -N 4 --export=MODE=prod TES_NORTHCanESM5forcingGEN.sub
```

### Commands and instructions

- Quick test (default case index 0, 1980-1999, time_steps=1):
```
sbatch /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_generation_canESM/TES_NORTHCanESM5forcingGEN.sub
```
- Quick test with overrides (single case, 20 years):
```
sbatch --export=MODE=test,TEST_CASE_INDEX=2,TEST_START_YEAR=1980,TEST_YEAR_COUNT=20 \
  /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_generation_canESM/TES_NORTHCanESM5forcingGEN.sub
```
- Production run (4 nodes, all timesteps):
```
sbatch -N 4 --export=MODE=prod \
  /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_generation_canESM/TES_NORTHCanESM5forcingGEN.sub
```
- Slurm output appears as `slurm-<jobid>.out` in the submission directory unless redirected by Slurm defaults.
- Case logs are written to this folder as `<EXPID>_forcinggen*.log.<timestamp>`.

### Output locations

- Test outputs: `/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_cases_data/<EXPID>/forcing_test_<start>_<end>/`
- Production outputs: `/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_cases_data/<EXPID>/forcing/`

Logs are written under this folder as:
`<EXPID>_forcinggen*.log.<timestamp>`.
