# Portable QEPP build
# Builds against Ubuntu 20.04 (GLIBC 2.31) so the binary runs on any
# modern Linux without requiring a matching GLIBC or libstdc++ version.
#
# Usage:
#   docker build -t qepp-builder .
#   docker run --rm -v "$PWD/build-docker:/out" qepp-builder \
#       cp /src/build/qepp /out/qepp

FROM ubuntu:20.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        libeigen3-dev \
        gnuplot \
        libfftw3-dev \
        libpng-dev \
        libjpeg-dev \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ── Build Matplot++ 1.2.0 from source ──────────────────────────────────────
RUN git clone --depth 1 --branch v1.2.0 \
        https://github.com/alandefreitas/matplotplusplus.git /tmp/matplot && \
    cmake -S /tmp/matplot -B /tmp/matplot/build \
          -DCMAKE_BUILD_TYPE=Release \
          -DMATPLOTPP_BUILD_EXAMPLES=OFF \
          -DMATPLOTPP_BUILD_TESTS=OFF && \
    cmake --build /tmp/matplot/build -j"$(nproc)" && \
    cmake --install /tmp/matplot/build && \
    rm -rf /tmp/matplot

# ── Build QEPP ─────────────────────────────────────────────────────────────
WORKDIR /src
COPY . .

RUN cmake -S . -B build \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_EXPORT_COMPILE_COMMANDS=ON && \
    cmake --build build -j"$(nproc)"

# ── Minimal output image containing only the binary ────────────────────────
FROM scratch AS export
COPY --from=builder /src/build/qepp /qepp
