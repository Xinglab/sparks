#!/bin/bash

function set_script_dir() {
  local ORIG_DIR="$(pwd)" || return 1

  local REL_SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})" || return 1
  cd "${REL_SCRIPT_DIR}" || return 1
  SCRIPT_DIR="$(pwd)" || return 1
  cd "${ORIG_DIR}" || return 1
}

function install_r_package() {
  Rscript install_package.R || return 1
}

function install_conda() {
  conda create --prefix "${CONDA_ENV_PREFIX}" || return 1
  conda activate "${CONDA_ENV_PREFIX}" || return 1
  conda install -c conda-forge -c bioconda --file \
        "${SCRIPT_DIR}/conda_requirements.txt" || return 1
  install_r_package || return 1
  conda deactivate || return 1
}

function write_config_paths() {
  echo "conda_wrapper: '${SCRIPT_DIR}/conda_wrapper'" \
       >> "${SCRIPT_DIR}/snakemake_config.yaml" || return 1
  echo "scripts_dir: '${SCRIPT_DIR}/aux_scripts'" \
       >> "${SCRIPT_DIR}/snakemake_config.yaml" || return 1
}

function install() {
  install_conda || return 1
  write_config_paths || return 1
}

function main() {
  set_script_dir || return 1
  source "${SCRIPT_DIR}/set_env_vars.sh" || return 1

  install || return 1
}

main "$@"
