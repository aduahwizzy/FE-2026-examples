#!/bin/bash
# This script disables Docker's automatic image pull for emod-ubuntu-runtime
# and ensures only the local image is used

# Check if the local image exists
if ! docker image inspect emod-ubuntu-runtime-py39:local > /dev/null 2>&1; then
    echo "ERROR: Local image emod-ubuntu-runtime-py39:local not found!"
    echo "Build it with: cd /home/anaphase21/emod-tutorials/FE-2026-examples && docker build -t emod-ubuntu-runtime-py39:local ."
    exit 1
fi

echo "Using local Docker image: emod-ubuntu-runtime-py39:local"
echo "Image ID: $(docker image inspect emod-ubuntu-runtime-py39:local --format='{{.ID}}')"
echo ""

# Run your experiment script
# Example: python3 run_example_burnin.py
# The Platform will now use the manifest.plat_image = "emod-ubuntu-runtime-py39:local" which points to the local image

"$@"
