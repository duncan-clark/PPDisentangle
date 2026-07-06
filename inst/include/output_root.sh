# Shared helper for bash drivers. Source after PKG_ROOT is set.
#
#   source "$PKG_ROOT/inst/include/output_root.sh"
#   OUTPUT_DIR="$(pp_disentangle_output_path oklahoma)"
#
# Override the default sibling folder with PPDISENTANGLE_OUTPUT_ROOT.

pp_disentangle_output_root() {
  if [ -n "${PPDISENTANGLE_OUTPUT_ROOT:-}" ]; then
    printf '%s\n' "$PPDISENTANGLE_OUTPUT_ROOT"
  else
    printf '%s\n' "$(dirname "$PKG_ROOT")/PPDisentangle-output"
  fi
}

pp_disentangle_output_path() {
  printf '%s/%s\n' "$(pp_disentangle_output_root)" "$1"
}
