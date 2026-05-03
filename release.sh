#!/usr/bin/env bash
# release.sh — build qepp in Release mode and optionally push source to GitHub
# The Release build always links spglib statically (via FetchContent if needed),
# so the resulting binary has no runtime dependency on libsymspg.so.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -b, --build          Build the release binary (default if no options given)
  -p, --push [MSG]     Commit any staged/unstaged changes and push to origin/main
                       MSG is optional; defaults to "Update source"
  -a, --all [MSG]      Build then push (equivalent to -b -p)
  -h, --help           Show this help

Examples:
  $0                          # build only
  $0 --push "Fix bug in dos"  # commit + push only
  $0 --all "Release v1.2"     # build then commit + push
EOF
}

do_build=0
do_push=0
commit_msg="Update source"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--build)
            do_build=1
            shift
            ;;
        -p|--push)
            do_push=1
            shift
            if [[ $# -gt 0 && "$1" != -* ]]; then
                commit_msg="$1"
                shift
            fi
            ;;
        -a|--all)
            do_build=1
            do_push=1
            shift
            if [[ $# -gt 0 && "$1" != -* ]]; then
                commit_msg="$1"
                shift
            fi
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

# Default to build if nothing specified
if [[ $do_build -eq 0 && $do_push -eq 0 ]]; then
    do_build=1
fi

# ── Build ──────────────────────────────────────────────────────────────────
if [[ $do_build -eq 1 ]]; then
    echo "==> Configuring release build..."
    cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

    echo "==> Building..."
    cmake --build "$BUILD_DIR" -j "$(nproc)"

    echo "==> Build complete: $BUILD_DIR/qepp"
fi

# ── Push ───────────────────────────────────────────────────────────────────
if [[ $do_push -eq 1 ]]; then
    cd "$SCRIPT_DIR"

    if git diff --quiet && git diff --cached --quiet; then
        echo "==> Nothing to commit."
    else
        echo "==> Staging all changes..."
        git add -A
        echo "==> Committing: \"$commit_msg\""
        git commit -m "$commit_msg"
    fi

    LOCAL=$(git rev-parse @)
    REMOTE=$(git rev-parse "@{u}" 2>/dev/null || echo "")
    if [[ -z "$REMOTE" || "$LOCAL" != "$REMOTE" ]]; then
        echo "==> Pushing to origin..."
        git push origin main
        echo "==> Push complete."
    else
        echo "==> Already up to date with origin; nothing to push."
    fi
fi
