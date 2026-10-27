#!/bin/bash
# Build the PHREEQC command-line program from the USGS source.
# ARTEMIS v1.1's published outputs were made with PHREEQC 3.7.3; later
# versions (3.8.x, 3.9.0) differ by up to ~0.01 pH after kinetic steps.
#
# Usage: scripts/build_phreeqc.sh [tag] [dest]   (defaults: v3.7.3, ./phreeqc-<tag>)
set -euo pipefail
tag=${1:-v3.7.3}
dest=${2:-phreeqc-$tag}
git clone --depth 1 --branch "$tag" https://github.com/usgs-coupled/phreeqc3 "$dest"
cmake -S "$dest" -B "$dest/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$dest/build" -j"$(nproc)"
echo "PHREEQC built: $dest/build/phreeqc  (export PHREEQC=$(cd "$dest/build" && pwd)/phreeqc)"
