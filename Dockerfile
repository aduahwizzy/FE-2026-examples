# Docker configuration for emod-malaria 5.0.2 with Python 3.9 and MPI support
# 
# The emod-malaria 5.0.2 binary is compiled against:
# - Python 3.9 (runtime libraries)
# - MPICH (MPI library with symlinks for libmpi.so.12 and libmpicxx.so.12)
#
# This Dockerfile extends the runtime image to add the required libraries.
# Build: docker build -t ghcr.io/emod-hub/emod-ubuntu-runtime:5.0.2-py39 .

FROM ghcr.io/emod-hub/emod-ubuntu-runtime:latest@sha256:91933ca254ac9c0dd49deb6ba9b48c59b312c2baae2970e547b6a6f5b896fbbd

# Install Python 3.9 runtime libraries required by emod-malaria 5.0.2 binary
RUN apt-get update && apt-get install -y \
    libpython3.9 \
    libpython3.9-stdlib \
    libpython3.9-minimal \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Fix MPI library symlinks: create symlinks for libmpi.so.12 and libmpicxx.so.12
# The emod-malaria binary looks for these but the system has libmpich* instead
RUN ln -sf /usr/lib/x86_64-linux-gnu/libmpich.so.12 /usr/lib/x86_64-linux-gnu/libmpi.so.12 && \
    ln -sf /usr/lib/x86_64-linux-gnu/libmpichcxx.so.12 /usr/lib/x86_64-linux-gnu/libmpicxx.so.12

# Verify installation
RUN echo "Verifying Python 3.9..." && \
    find /usr/lib -name "libpython3.9.so*" > /dev/null && echo "✓ Python 3.9 libraries installed" && \
    echo "Verifying MPI..." && \
    ls -la /usr/lib/x86_64-linux-gnu/libmpi.so.12 && \
    ls -la /usr/lib/x86_64-linux-gnu/libmpicxx.so.12 && \
    echo "✓ MPI library symlinks created"
