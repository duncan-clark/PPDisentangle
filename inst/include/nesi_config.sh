# NeSI SSH/rsync configuration for local launchers under inst/nesi/.
#
# Copy inst/nesi/config.example.env to inst/nesi/nesi.env (gitignored) or set
# PPDISENTANGLE_NESI_CONFIG to your own env file path.

pp_nesi_load_config() {
  local pkg_root="$1"
  local candidates=()
  if [ -n "${PPDISENTANGLE_NESI_CONFIG:-}" ]; then
    candidates+=("$PPDISENTANGLE_NESI_CONFIG")
  fi
  candidates+=(
    "$HOME/.config/ppdisentangle/nesi.env"
    "$pkg_root/inst/nesi/nesi.env"
  )
  local f
  for f in "${candidates[@]}"; do
    if [ -f "$f" ]; then
      # shellcheck disable=SC1090
      source "$f"
      echo "Loaded NeSI config: $f"
      return 0
    fi
  done
  return 0
}

pp_nesi_apply_defaults() {
  local pkg_root="$1"
  : "${PP_NESI_USER:=${USER:-$(id -un)}}"
  : "${PP_NESI_HOST:=mahuika.nesi.org.nz}"
  : "${PP_NESI_SSH:=${PP_NESI_USER}@${PP_NESI_HOST}}"
  : "${PP_NESI_REMOTE_PKG:=/nesi/home/users/${PP_NESI_USER}/PPDisentangle}"
  : "${PP_NESI_REMOTE_OUTPUT:=/nesi/project/uoo04008/PPDisentangle-output}"

  export PP_NESI_USER PP_NESI_HOST PP_NESI_SSH
  export PP_NESI_REMOTE_PKG PP_NESI_REMOTE_OUTPUT

  # Local output root (sibling to this checkout by default).
  # shellcheck source=inst/include/output_root.sh
  source "$pkg_root/inst/include/output_root.sh"
  PKG_ROOT="$pkg_root"
  export PP_NESI_LOCAL_OUTPUT="$(pp_disentangle_output_root)"
  export PKG_ROOT
}

pp_nesi_init_config() {
  local pkg_root="$1"
  pp_nesi_load_config "$pkg_root"
  pp_nesi_apply_defaults "$pkg_root"
}

pp_nesi_print_config() {
  cat <<EOF
NeSI configuration:
  PP_NESI_SSH              = ${PP_NESI_SSH}
  PP_NESI_REMOTE_PKG       = ${PP_NESI_REMOTE_PKG}
  PP_NESI_REMOTE_OUTPUT    = ${PP_NESI_REMOTE_OUTPUT}
  PP_NESI_LOCAL_OUTPUT     = ${PP_NESI_LOCAL_OUTPUT}
  PP_GIT_BRANCH            = ${PP_GIT_BRANCH:-<remote/current>}
  PP_NESI_PUSH             = ${PP_NESI_PUSH:-0}
EOF
}
