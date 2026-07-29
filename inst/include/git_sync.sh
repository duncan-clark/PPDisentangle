# Git sync helpers for NeSI launchers and local submit wrappers.
#
#   source "$PKG_ROOT/inst/include/git_sync.sh"
#   pp_git_sync_repo "$PKG_ROOT"

pp_git_current_branch() {
  git -C "$1" rev-parse --abbrev-ref HEAD 2>/dev/null || true
}

pp_git_sync_repo() {
  local repo_root="$1"
  if ! git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Not a git repo; skipping sync."
    return 0
  fi
  if [ "${PP_GIT_SYNC:-1}" = "0" ]; then
    echo "Skipping git sync (PP_GIT_SYNC=0)."
    return 0
  fi

  local branch="${PP_GIT_BRANCH:-}"
  if [ -z "$branch" ]; then
    branch="$(pp_git_current_branch "$repo_root")"
  fi
  if [ -z "$branch" ] || [ "$branch" = "HEAD" ]; then
    echo "Detached HEAD and PP_GIT_BRANCH unset; skipping git sync."
    return 0
  fi

  echo "Git sync: fetch/pull origin/$branch in $repo_root"
  git -C "$repo_root" fetch origin "$branch" 2>/dev/null || git -C "$repo_root" fetch origin 2>/dev/null || true
  if git -C "$repo_root" pull --ff-only origin "$branch" 2>/dev/null; then
    return 0
  fi
  if git -C "$repo_root" pull --ff-only 2>/dev/null; then
    return 0
  fi
  echo "Warning: git pull failed; continuing with the current checkout."
}

pp_git_push_branch() {
  local repo_root="$1"
  if [ "${PP_NESI_PUSH:-0}" != "1" ]; then
    return 0
  fi
  if ! git -C "$repo_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Not a git repo; skipping push."
    return 0
  fi

  local branch="${PP_GIT_BRANCH:-}"
  if [ -z "$branch" ]; then
    branch="$(pp_git_current_branch "$repo_root")"
  fi
  if [ -z "$branch" ] || [ "$branch" = "HEAD" ]; then
    echo "Cannot push: detached HEAD. Set PP_GIT_BRANCH."
    return 1
  fi

  echo "Git push: origin $branch from $repo_root"
  git -C "$repo_root" push origin "$branch"
}
