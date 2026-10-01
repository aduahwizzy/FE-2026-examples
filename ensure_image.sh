#!/bin/bash
# Ensure the Docker image is correctly set up before running EMOD experiments
# This prevents Docker Desktop from reverting the 'latest' tag to the old image

set -e

IMAGE_TAG="ghcr.io/emod-hub/emod-ubuntu-runtime:latest"
FIXED_IMAGE_TAG="ghcr.io/emod-hub/emod-ubuntu-runtime:5.0.2-py39"

echo "Checking Docker image setup for emodpy-malaria 5.0.2..."

# Check if the fixed image exists
if ! docker image inspect "$FIXED_IMAGE_TAG" > /dev/null 2>&1; then
    echo "ERROR: Fixed image $FIXED_IMAGE_TAG not found!"
    echo "Rebuilding..."
    docker build -t "$FIXED_IMAGE_TAG" /home/anaphase21/emod-tutorials/FE-2026-examples/
fi

# Check if latest tag points to the fixed image
LATEST_ID=$(docker inspect "$IMAGE_TAG" --format='{{.ID}}' 2>/dev/null || echo "not-found")
FIXED_ID=$(docker inspect "$FIXED_IMAGE_TAG" --format='{{.ID}}')

if [ "$LATEST_ID" != "$FIXED_ID" ]; then
    echo "Retagging $IMAGE_TAG to point to $FIXED_IMAGE_TAG..."
    docker rmi -f "$IMAGE_TAG" 2>/dev/null || true
    docker tag "$FIXED_IMAGE_TAG" "$IMAGE_TAG"
fi

# Verify Python 3.9 is present
echo "Verifying Python 3.9 libraries..."
docker run --rm "$IMAGE_TAG" bash -c "ls /usr/lib/x86_64-linux-gnu/libpython3.9.so.1.0 > /dev/null && echo '✓ Python 3.9 ready'" || {
    echo "ERROR: Python 3.9 not found in image!"
    exit 1
}

echo "✓ Docker image is ready for emodpy-malaria 5.0.2"
echo ""
echo "You can now run your experiment:"
echo "  python3 run_example_pickup_CM.py"
