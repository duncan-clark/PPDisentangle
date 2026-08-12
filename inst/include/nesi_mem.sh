# Memory helpers for NeSI launchers / run scripts.
# Bootstrap outer workers are mostly single-threaded and memory-heavy
# (Oklahoma bootstrap-only ~0.85 GB RSS/worker; full-mode spikes ~1.9 GB).

pp_mem_to_gb() {
  local raw="${1:-0}"
  raw="$(printf '%s' "$raw" | tr '[:lower:]' '[:upper:]' | tr -d ' ')"
  local num="${raw%%[A-Z]*}"
  local unit="${raw#"$num"}"
  if [ -z "$num" ]; then
    echo 0
    return 0
  fi
  case "$unit" in
    M|MB) echo $(( (num + 1023) / 1024 )) ;;
    T|TB) echo $(( num * 1024 )) ;;
    G|GB|"") echo "$num" ;;
    *) echo "$num" ;;
  esac
}

# Cap outer PSOCK workers so cores * worker_gb + reserve fits in PP_MEM.
# Usage: pp_boot_outer_from_mem CORES MEM [WORKER_GB] [RESERVE_GB]
pp_boot_outer_from_mem() {
  local cores="${1:-1}"
  local mem="${2:-64G}"
  local worker_gb="${3:-2}"
  local reserve_gb="${4:-8}"
  local mem_gb ram_cap
  mem_gb="$(pp_mem_to_gb "$mem")"
  if [ "$mem_gb" -lt 1 ]; then mem_gb=1; fi
  if [ "$worker_gb" -lt 1 ]; then worker_gb=1; fi
  if [ "$reserve_gb" -lt 0 ]; then reserve_gb=0; fi
  ram_cap=$(( (mem_gb - reserve_gb) / worker_gb ))
  if [ "$ram_cap" -lt 1 ]; then ram_cap=1; fi
  if [ "$cores" -gt 0 ] && [ "$ram_cap" -gt "$cores" ]; then
    ram_cap="$cores"
  fi
  echo "$ram_cap"
}
