#!/bin/bash

EXPID="CanESM5_1pctCO2-bgc_r1i1p1f1_DBCCA_Daymet_TESSFA1"
SRC="/gpfs/wolf2/cades/cli185/proj-shared/wangd/canESM5/${EXPID}"
SUBSET="/gpfs/wolf2/cades/cli185/proj-shared/wangd/canESM5/${EXPID}_subset_TPHWL3Hrly_Solar3Hrly"

mkdir -p "${SUBSET}"
ln -sfn "${SRC}/TPHWL3Hrly" "${SUBSET}/TPHWL3Hrly"
ln -sfn "${SRC}/Solar3Hrly" "${SUBSET}/Solar3Hrly"

srun -N 1 -n 1 --time=12:00:00 \
  /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/python_test_env/conda_envs/testvenv/bin/python \
  /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_inputGEN/TES_NORTH/data_generation_canESM/TES_forcingGEN_NORTHCanESM5_1pctCO2-bgc_r1i1p1f1_DBCCA_Daymet_TESSFA1.py \
  "${SUBSET}" \
  "/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_cases_data/${EXPID}/forcing" \
  -1
