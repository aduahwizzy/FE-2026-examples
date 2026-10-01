# EMOD-Malaria Docker Fix for Python 3.9

## Problem
When running emod-malaria 5.0.2 experiments, you get:
```
error while loading shared libraries: libpython3.9.so.1.0: cannot open shared object file: No such file or directory
```

This occurs because the emod-malaria binary is compiled against Python 3.9, but the official Docker image only includes Python 3.13.

## Solution
A custom Docker image has been built locally with Python 3.9 support added. All run scripts have been configured to use this image.

### Image Details
- **Local image name**: `emod-python39:local`
- **Base image**: `ghcr.io/emod-hub/emod-ubuntu-runtime:latest`
- **Additional packages**: libpython3.9, libpython3.9-stdlib, libpython3.9-minimal

### How It Works
1. The `Dockerfile` in this directory builds the custom image
2. `manifest.py` is configured to use `emod-python39:local` as the Docker image
3. All `run_example_*.py` scripts pass `plat_image=manifest.plat_image` to the Platform
4. When you run experiments, they use the local image instead of pulling from the registry

### Verification
Your local image is ready. Run any experiment script:
```bash
python3 run_example_burnin.py
```

The framework will use the local `emod-python39:local` image, which has Python 3.9 support built in.

### If You Need to Rebuild
If you ever need to rebuild the image:
```bash
cd /home/anaphase21/emod-tutorials/FE-2026-examples
docker build -t emod-python39:local .
```

### Files Modified
- `Dockerfile` - New: builds the custom image with Python 3.9
- `manifest.py` - Updated: changed `plat_image` to use `emod-python39:local`
- `run_example*.py` - Updated: added `plat_image=manifest.plat_image` to Platform initialization in all 6 files:
  - run_example.py
  - run_example_burnin.py
  - run_example_calibration.py
  - run_example_outputs.py
  - run_example_pickup.py
  - run_example_sweeps.py
