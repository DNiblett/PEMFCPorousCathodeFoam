#!/usr/bin/env bash
# Build the tutorial applications into $FOAM_USER_APPBIN.
set -euo pipefail

if ! command -v wmake >/dev/null 2>&1; then
    echo "Source the OpenFOAM-v2512 etc/bashrc before building." >&2
    exit 1
fi

root="$(cd "$(dirname "$0")/.." && pwd)"
for application in pemfcPorousFlowFoam pemfcPorousCathodeFoam; do
    printf 'Building %s\n' "$application"
    (cd "$root/applications/$application" && wmake)
done
