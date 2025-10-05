#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import stat
import sys
from datetime import datetime
from pathlib import Path


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def make_executable(path: Path) -> None:
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_text_file(path: Path, content: str) -> None:
    path.write_text(content)


def copy_if_missing(src: Path, dst: Path) -> None:
    if not dst.exists():
        shutil.copy2(src, dst)


def resolve_required_file(path_str: str, description: str) -> Path:
    p = Path(path_str).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"Required {description} not found: {p}")
    return p


def render_run_domain_surfdata_sh(cfg: dict, scripts_dir: Path, exp_root: Path) -> str:
    expid = cfg["expid"]
    aoi_dir = cfg["aoi_points"]["dir"].rstrip("/")
    aoi_file = cfg["aoi_points"]["file"]
    base_domain_file = cfg["source"]["base_domain_file"]
    surf_dir = cfg["source"]["surfdata_dir"].rstrip("/")
    surf_file = cfg["source"]["surfdata_file"]

    scheduler = cfg.get("scheduler", {})
    module_line = scheduler.get("module", "")
    conda_activate = scheduler.get("conda_activate", "")

    lines = []
    lines.append("#!/bin/bash")
    lines.append("set -euo pipefail")
    if module_line:
        lines.append(f"module load {module_line}")
    if conda_activate:
        # Use new 'conda activate' if string seems like env name, otherwise 'source activate' for explicit path
        if os.path.sep in conda_activate:
            lines.append(f"source activate {conda_activate}")
        else:
            lines.append(f"conda activate {conda_activate}")
    lines.append("")
    lines.append(f"cd '{scripts_dir.as_posix()}'")
    lines.append("# Source exported environment if present")
    lines.append("if [ -f ./export_env.sh ]; then . ./export_env.sh; fi")
    lines.append("")
    lines.append("date_string=$(date +'%y%m%d-%H%M')")
    lines.append(f": \"${{EXPID:={expid}}}\"")
    lines.append(f": \"${{EXP_ROOT:={exp_root.as_posix()}}}\"")
    lines.append(f": \"${{AOI_POINTS_DIR:={aoi_dir}}}\"")
    lines.append(f": \"${{AOI_POINTS_FILE:={aoi_file}}}\"")
    lines.append(f": \"${{BASE_DOMAIN_FILE:={base_domain_file}}}\"")
    lines.append(f": \"${{SURFDATA_DIR:={surf_dir}}}\"")
    lines.append(f": \"${{SURFDATA_FILE:={surf_file}}}\"")
    lines.append("DOM_SURF_DIR=\"${EXP_ROOT}/domain_surfdata\"")
    lines.append("mkdir -p \"${DOM_SURF_DIR}\"")
    lines.append("")
    lines.append("# Ensure domain source file is available with expected local name")
    lines.append("if [ ! -e domain.lnd.TES_SE.4km.1d.nc ]; then")
    lines.append("  ln -sf \"${BASE_DOMAIN_FILE}\" domain.lnd.TES_SE.4km.1d.nc")
    lines.append("fi")
    lines.append("")
    lines.append("echo \"[1/2] Generating AOI domain...\"")
    lines.append("python3 TES_AOI_domainGEN.py \"${AOI_POINTS_DIR}\" \"${DOM_SURF_DIR}\" \"${AOI_POINTS_FILE}\" 2>&1 | tee \"${DOM_SURF_DIR}/${EXPID}_domaingen.log.${date_string}\"")
    lines.append("")
    lines.append("echo \"Resolving latest AOI domain file...\"")
    lines.append("AOI_PREFIX=\"$(echo \"${AOI_POINTS_FILE}\" | awk -F'_' '{print $1}')\"")
    lines.append("AOI_DOMAIN=$(ls -1 ${DOM_SURF_DIR}/${AOI_PREFIX}_domain.lnd.TES_SE.4km.1d.c*.nc 2>/dev/null | sort | tail -n1)")
    lines.append("if [ -z \"${AOI_DOMAIN}\" ]; then echo 'ERROR: AOI domain file not found'; exit 2; fi")
    lines.append("")
    lines.append("echo \"[2/2] Generating AOI surfdata...\"")
    lines.append("python3 TES_AOI_surfdataGEN.py \"${SURFDATA_DIR}\" \"${SURFDATA_FILE}\" \"${DOM_SURF_DIR}\" \"${DOM_SURF_DIR}/\" \"$(basename \"${AOI_DOMAIN}\")\" 2>&1 | tee \"${DOM_SURF_DIR}/${EXPID}_surfdargen.log.${date_string}\"")
    lines.append("")
    lines.append("echo 'Domain and surfdata generation complete.'")
    return "\n".join(lines) + "\n"


def render_run_forcing_sbatch(cfg: dict, scripts_dir: Path, exp_root: Path) -> str:
    expid = cfg["expid"]
    forcing_dir = cfg["source"]["forcing_dir"].rstrip("/")
    scheduler = cfg.get("scheduler", {})
    account = scheduler.get("account", "")
    partition = scheduler.get("partition", "batch")
    nodes = scheduler.get("nodes", 1)
    time_limit = scheduler.get("time", "24:00:00")
    mem = scheduler.get("mem", "128GB")
    module_line = scheduler.get("module", "")
    conda_activate = scheduler.get("conda_activate", "")

    lines = []
    lines.append("#!/bin/bash")
    if account:
        lines.append(f"#SBATCH -A {account}")
    lines.append(f"#SBATCH -J TES_{expid}_forcingGEN")
    lines.append(f"#SBATCH -p {partition}")
    lines.append(f"#SBATCH -N {nodes}")
    lines.append(f"#SBATCH -t {time_limit}")
    lines.append(f"#SBATCH --mem={mem}")
    lines.append("")
    lines.append("set -euo pipefail")
    if module_line:
        lines.append(f"module load {module_line}")
    if conda_activate:
        if os.path.sep in conda_activate:
            lines.append(f"source activate {conda_activate}")
        else:
            lines.append(f"conda activate {conda_activate}")
    lines.append("")
    lines.append(f"cd '{scripts_dir.as_posix()}'")
    lines.append("# Source exported environment if present")
    lines.append("if [ -f ./export_env.sh ]; then . ./export_env.sh; fi")
    lines.append("")
    lines.append("date_string=$(date +'%y%m%d-%H%M')")
    lines.append(f": \"${{EXPID:={expid}}}\"")
    lines.append(f": \"${{EXP_ROOT:={exp_root.as_posix()}}}\"")
    lines.append(f": \"${{FORCING_DIR:={forcing_dir}}}\"")
    lines.append("OUT_DIR=\"${EXP_ROOT}/forcing\"")
    lines.append("mkdir -p \"${OUT_DIR}\"")
    lines.append("")
    lines.append("# Locate the latest AOI domain to derive gridIDs")
    lines.append("AOI_FILE_PATH=\"${EXP_ROOT}/domain_surfdata\"")
    lines.append("AOI_POINTS_FILE=$(ls -1 ${AOI_FILE_PATH}/*_domain.lnd.TES_SE.4km.1d.c*.nc 2>/dev/null | sort | tail -n1 | xargs -r basename)")
    lines.append("if [ -z \"${AOI_POINTS_FILE}\" ]; then echo 'ERROR: AOI domain file not found'; exit 2; fi")
    lines.append("")
    lines.append("# If not running inside a Slurm allocation, run with srun using env SCHED_* overrides")
    lines.append("if [ -z \"${SLURM_JOB_ID:-}\" ]; then")
    lines.append("  ACCOUNT=\"${SCHED_ACCOUNT:-" + account + "}\"")
    lines.append("  PARTITION=\"${SCHED_PARTITION:-" + partition + "}\"")
    lines.append("  NODES=\"${SCHED_NODES:-" + str(nodes) + "}\"")
    lines.append("  TIME=\"${SCHED_TIME:-" + time_limit + "}\"")
    lines.append("  MEM=\"${SCHED_MEM:-" + mem + "}\"")
    lines.append("  echo \"srun -A '${ACCOUNT}' -p '${PARTITION}' -N '${NODES}' -t '${TIME}' --mem='${MEM}' python3 TES_AOI_forcingGEN.py '${FORCING_DIR}' '${OUT_DIR}' '${AOI_FILE_PATH}/' '${AOI_POINTS_FILE}'\" | tee \"${OUT_DIR}/${EXPID}_forcinggen.cmd.${date_string}\"")
    lines.append("  exec srun -A \"${ACCOUNT}\" -p \"${PARTITION}\" -N \"${NODES}\" -t \"${TIME}\" --mem=\"${MEM}\" python3 TES_AOI_forcingGEN.py \"${FORCING_DIR}\" \"${OUT_DIR}\" \"${AOI_FILE_PATH}/\" \"${AOI_POINTS_FILE}\" 2>&1 | tee \"${OUT_DIR}/${EXPID}_forcinggen.log.${date_string}\"")
    lines.append("fi")
    lines.append("")
    lines.append("# Running under Slurm allocation")
    lines.append("echo \"python3 TES_AOI_forcingGEN.py '${FORCING_DIR}' '${OUT_DIR}' '${AOI_FILE_PATH}/' '${AOI_POINTS_FILE}'\" | tee \"${OUT_DIR}/${EXPID}_forcinggen.cmd.${date_string}\"")
    lines.append("srun -n 1 python3 TES_AOI_forcingGEN.py \"${FORCING_DIR}\" \"${OUT_DIR}\" \"${AOI_FILE_PATH}/\" \"${AOI_POINTS_FILE}\" 2>&1 | tee \"${OUT_DIR}/${EXPID}_forcinggen.log.${date_string}\"")
    return "\n".join(lines) + "\n"


def render_create_links_sh() -> str:
    lines = []
    lines.append("#!/bin/bash")
    lines.append("cd \"$(dirname \"$0\")\"")
    lines.append("python3 forcing_domain_link_creation.py")
    lines.append("echo 'Soft links created under ../atm_forcing.datm7.km.1d' ")
    return "\n".join(lines) + "\n"


def render_export_env_sh(cfg: dict, exp_root: Path) -> str:
    expid = cfg["expid"]
    aoi_dir = cfg["aoi_points"]["dir"].rstrip("/")
    aoi_file = cfg["aoi_points"]["file"]
    base_domain_file = cfg["source"]["base_domain_file"]
    surf_dir = cfg["source"]["surfdata_dir"].rstrip("/")
    surf_file = cfg["source"]["surfdata_file"]
    forcing_dir = cfg["source"]["forcing_dir"].rstrip("/")

    scheduler = cfg.get("scheduler", {})

    lines = []
    lines.append("#!/bin/bash")
    lines.append("# Auto-generated by aoi_prepare_experiment.py from your config")
    lines.append(f"export EXPID=\"{expid}\"")
    lines.append(f"export EXP_ROOT=\"{exp_root.as_posix()}\"")
    lines.append(f"export AOI_POINTS_DIR=\"{aoi_dir}\"")
    lines.append(f"export AOI_POINTS_FILE=\"{aoi_file}\"")
    lines.append(f"export BASE_DOMAIN_FILE=\"{base_domain_file}\"")
    lines.append(f"export SURFDATA_DIR=\"{surf_dir}\"")
    lines.append(f"export SURFDATA_FILE=\"{surf_file}\"")
    lines.append(f"export FORCING_DIR=\"{forcing_dir}\"")
    lines.append("")
    # Optional scheduler exports for external use
    if scheduler:
        for key, value in scheduler.items():
            var = f"SCHED_{key.upper()}"
            lines.append(f"export {var}=\"{value}\"")
    return "\n".join(lines) + "\n"


def render_create_uelm_adspin_sh(cfg: dict, exp_root: Path) -> str:
    expid = cfg["expid"]
    e3sm = cfg.get("e3sm", {})
    din_root = e3sm.get("din_root", "//gpfs/wolf2/cades/cli185/world-shared/e3sm")
    src_root = e3sm.get("src_root", "$(git rev-parse --show-toplevel)")
    mach = e3sm.get("mach", "cades-baseline")
    compiler = e3sm.get("compiler", "gnu")
    mpilib = e3sm.get("mpilib", "openmpi")
    compset = e3sm.get("compset", "I1850uELMTESCNPRDCTCBC")

    lines = []
    lines.append("#!/bin/bash")
    lines.append("set -e")
    lines.append("")
    lines.append(f": \"${{EXPID:={expid}}}\"")
    lines.append(f": \"${{EXP_ROOT:={exp_root.as_posix()}}}\"")
    lines.append(f": \"${{E3SM_DIN:={din_root}}}\"")
    lines.append(f": \"${{E3SM_SRCROOT:={src_root}}}\"")
    lines.append("")
    lines.append("CASE_DATA=\"${EXP_ROOT}\"")
    lines.append("CASEDIR=\"${E3SM_SRCROOT}/e3sm_cases/uELM_${EXPID}_I1850uELMCNPRDCTCBC\"")
    lines.append("")
    # Discover latest domain/surfdata filenames
    lines.append("DOMAIN_FILE=$(ls -1 ${CASE_DATA}/domain_surfdata/${EXPID}_domain.lnd.TES_SE.4km.1d.c*.nc 2>/dev/null | sort | tail -n1 | xargs -r basename)")
    lines.append("if [ -z \"${DOMAIN_FILE}\" ]; then echo 'ERROR: Domain file not found'; exit 2; fi")
    lines.append("SURFDATA_FILE=$(ls -1 ${CASE_DATA}/domain_surfdata/${EXPID}_surfdata.TES_SE.4km.1d.NLCD.c*.nc 2>/dev/null | sort | tail -n1 | xargs -r basename)")
    lines.append("if [ -z \"${SURFDATA_FILE}\" ]; then echo 'ERROR: Surfdata file not found'; exit 2; fi")
    lines.append("")
    lines.append("rm -rf \"${CASEDIR}\"")
    lines.append("")
    lines.append(f"${{E3SM_SRCROOT}}/cime/scripts/create_newcase --case \"${{CASEDIR}}\" --mach {mach} --compiler {compiler} --mpilib {mpilib} --compset {compset} --res ELM_USRDAT  --handle-preexisting-dirs r --srcroot \"${{E3SM_SRCROOT}}\"")
    lines.append("")
    lines.append("cd \"${CASEDIR}\"")
    lines.append("")
    lines.append("./xmlchange PIO_TYPENAME=\"pnetcdf\"")
    lines.append("./xmlchange PIO_NETCDF_FORMAT=\"64bit_data\"")
    lines.append("./xmlchange DIN_LOC_ROOT=\"${E3SM_DIN}\"")
    lines.append("./xmlchange DIN_LOC_ROOT_CLMFORC=\"${CASE_DATA}\"")
    lines.append("./xmlchange CIME_OUTPUT_ROOT=\"${E3SM_SRCROOT}/e3sm_runs/\"")
    lines.append("./xmlchange ELM_FORCE_COLDSTART=on")
    lines.append("./xmlchange DATM_MODE=\"uELM_TES\"")
    lines.append("./xmlchange DATM_CLMNCEP_YR_START=\"1980\"")
    lines.append("./xmlchange DATM_CLMNCEP_YR_END=\"1999\"")


    lines.append("./xmlchange ATM_NCPL=\"24\"")
    lines.append("./xmlchange STOP_N=\"400\"")
    lines.append("./xmlchange STOP_OPTION=\"nyears\"")
    lines.append("./xmlchange NTASKS=\"1\"")
    lines.append("./xmlchange NTASKS_PER_INST=\"1\"")
    lines.append("./xmlchange MAX_MPITASKS_PER_NODE=\"128\"")
    lines.append("./xmlchange ATM_DOMAIN_PATH=\"${CASE_DATA}/domain_surfdata/\"")
    lines.append("./xmlchange ATM_DOMAIN_FILE=\"${DOMAIN_FILE}\"")
    lines.append("./xmlchange LND_DOMAIN_PATH=\"${CASE_DATA}/domain_surfdata/\"")
    lines.append("./xmlchange LND_DOMAIN_FILE=\"${DOMAIN_FILE}\"")
    lines.append("./xmlchange JOB_WALLCLOCK_TIME=\"2:00:00\"")
    lines.append("./xmlchange USER_REQUESTED_WALLTIME=\"2:00:00\"")
    lines.append("")

    lines.append("./xmlchange ELM_FORCE_COLDSTART=\"on\"")
    lines.append("./xmlchange CONTINUE_RUN=\"FALSE\"")
    lines.append("./xmlchange ELM_ACCELERATED_SPINUP=\"on\"")
    lines.append("./xmlchange  --append ELM_BLDNML_OPTS=\"-bgc_spinup on\"")

    lines.append("cat >> user_nl_elm <<'EOF'")
    lines.append("fsurdat = '${CASE_DATA}/domain_surfdata/${SURFDATA_FILE}'")
    lines.append("spinup_state = 1")
    lines.append("suplphos = 'ALL'")
    lines.append("hist_nhtfrq=-175200")
    lines.append("hist_mfilt=1")
    lines.append("nyears_ad_carbon_only = 25")
    lines.append("spinup_mortality_factor = 10")
    lines.append("EOF")

    lines.append("")
    lines.append("./case.setup --reset")
    lines.append("./case.setup")
    lines.append("")
    lines.append("./case.build --clean-all")
    lines.append("./case.build")
    lines.append("")
    # lines.append("./case.submit")  # left to the user
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a TES AOI experiment directory and wrappers.")
    parser.add_argument("--config", required=True, help="Path to JSON config file")
    parser.add_argument("--run-domain-surfdata", action="store_true", help="Run domain and surfdata generation now")
    parser.add_argument("--submit-forcing", action="store_true", help="Submit forcing generation job now (sbatch)")
    args = parser.parse_args()

    scripts_root = Path(__file__).resolve().parent

    cfg_path = resolve_required_file(args.config, "config JSON")
    cfg = json.loads(cfg_path.read_text())

    expid = cfg["expid"]
    exp_root = Path(cfg["experiment_root"]).expanduser()
    ensure_dir(exp_root)

    # Create expected subdirs
    domain_surf_dir = exp_root / "domain_surfdata"
    forcing_dir_out = exp_root / "forcing"
    user_scripts_dir = exp_root / "scripts"
    for d in (domain_surf_dir, forcing_dir_out, user_scripts_dir):
        ensure_dir(d)

    # Copy core generator scripts to user scripts dir
    core_scripts = [
        "TES_AOI_domainGEN.py",
        "TES_AOI_surfdataGEN.py",
        "TES_AOI_forcingGEN.py",
        "forcing_domain_link_creation.py",
        "forcinglink_creation.py",
        "check_nc_compression.py",
    ]
    for name in core_scripts:
        src = scripts_root / name
        if src.exists():
            copy_if_missing(src, user_scripts_dir / name)

    # Also copy original templates for reference if they exist
    optional_templates = ["domain_creation.sh", "forcing_creation.sbatch"]
    for name in optional_templates:
        src = scripts_root / name
        if src.exists():
            copy_if_missing(src, user_scripts_dir / (name + ".orig"))

    # Generate customized wrappers
    run_domain_surfdata = render_run_domain_surfdata_sh(cfg, user_scripts_dir, exp_root)
    write_text_file(user_scripts_dir / "run_domain_surfdata.sh", run_domain_surfdata)
    make_executable(user_scripts_dir / "run_domain_surfdata.sh")

    run_forcing_sbatch = render_run_forcing_sbatch(cfg, user_scripts_dir, exp_root)
    write_text_file(user_scripts_dir / "run_forcing.sbatch", run_forcing_sbatch)
    make_executable(user_scripts_dir / "run_forcing.sbatch")

    create_links_sh = render_create_links_sh()
    write_text_file(user_scripts_dir / "create_links.sh", create_links_sh)
    make_executable(user_scripts_dir / "create_links.sh")

    export_env_sh = render_export_env_sh(cfg, exp_root)
    write_text_file(user_scripts_dir / "export_env.sh", export_env_sh)
    make_executable(user_scripts_dir / "export_env.sh")

    # Add uELM adspin creator
    create_uelm_adspin = render_create_uelm_adspin_sh(cfg, exp_root)
    write_text_file(user_scripts_dir / "create_uELM_adspin.sh", create_uelm_adspin)
    make_executable(user_scripts_dir / "create_uELM_adspin.sh")

    # Optionally execute steps now
    if args.run_domain_surfdata:
        os.system(f"bash '{(user_scripts_dir / 'run_domain_surfdata.sh').as_posix()}'")

    if args.submit_forcing:
        os.system(f"sbatch '{(user_scripts_dir / 'run_forcing.sbatch').as_posix()}'")

    print("Prepared experiment at:", exp_root)
    print("- domain_surfdata:", domain_surf_dir)
    print("- forcing:", forcing_dir_out)
    print("- scripts:", user_scripts_dir)


if __name__ == "__main__":
    main()


