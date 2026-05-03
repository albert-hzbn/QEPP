#!/usr/bin/env bash
# release.sh — build qepp in Release mode and optionally push source to GitHub
# The Release build always links spglib and the C++ runtime statically.
# Use --docker to produce a binary that runs on systems with GLIBC >= 2.31.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
  -b, --build          Build the release binary (default if no options given)
      --docker         Build inside Ubuntu 20.04 container (GLIBC 2.31) for a
                       binary that runs on any Linux >= Ubuntu 20.04 / RHEL 8
  -p, --push [MSG]     Commit any staged/unstaged changes and push to origin/main
                       MSG is optional; defaults to "Update source"
  -a, --all [MSG]      Build then push (equivalent to -b -p)
  -h, --help           Show this help

Examples:
  $0                            # build only (native)
  $0 --docker                   # portable build via Docker (GLIBC 2.31)
  $0 --docker --all "v1.2"      # portable build then commit + push
  $0 --push "Fix bug in dos"    # commit + push only
  $0 --all "Release v1.2"       # native build then commit + push
EOF
}

do_build=0
do_push=0
do_docker=0
commit_msg="Update source"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--build)
            do_build=1
            shift
            ;;
        --docker)
            do_build=1
            do_docker=1
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
    if [[ $do_docker -eq 1 ]]; then
        # Build inside Ubuntu 20.04 container to target GLIBC 2.31
        if ! command -v docker &>/dev/null; then
            echo "ERROR: docker not found. Install Docker and retry." >&2
            exit 1
        fi
        DOCKER_OUT="$SCRIPT_DIR/build-docker"
        mkdir -p "$DOCKER_OUT"
        echo "==> Building portable binary inside Ubuntu 20.04 container..."
        docker build --target builder -t qepp-builder "$SCRIPT_DIR"
        docker run --rm -v "$DOCKER_OUT:/out" qepp-builder \
            cp /src/build/qepp /out/qepp
        echo "==> Portable build complete: $DOCKER_OUT/qepp"
    else
        echo "==> Configuring release build..."
        cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" \
              -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

        echo "==> Building..."
        cmake --build "$BUILD_DIR" -j "$(nproc)"

        echo "==> Build complete: $BUILD_DIR/qepp"
    fi
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
