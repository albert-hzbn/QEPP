# Installing QEPP

## Requirements

| Dependency | Minimum Version | Notes |
|------------|----------------|-------|
| C++ compiler | GCC 9 / Clang 10 | C++17 required |
| CMake | 3.16 | |
| [Eigen3](https://eigen.tuxfamily.org/) | 3.3 | Header-only; used for lattice / k-vector math |
| [Matplot++](https://alandefreitas.github.io/matplotplusplus/) | 1.1 | Plotting backend |
| [gnuplot](http://www.gnuplot.info/) | 5.4 | Runtime requirement for Matplot++ PNG output |
| [Quantum ESPRESSO](https://www.quantum-espresso.org/) | 7.x | `pw.x`, `dos.x`, `bands.x` must be on `PATH` |

---

## 1 — Install System Dependencies

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y build-essential cmake libeigen3-dev gnuplot
```

### Fedora / RHEL

```bash
sudo dnf install -y gcc-c++ cmake eigen3-devel gnuplot
```

### macOS (Homebrew)

```bash
brew install cmake eigen gnuplot
```

---

## 2 — Install Matplot++

Matplot++ is not available in most distro repositories; build it from source:

```bash
git clone https://github.com/alandefreitas/matplotplusplus.git
cd matplotplusplus
cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DMATPLOTPP_BUILD_EXAMPLES=OFF \
      -DMATPLOTPP_BUILD_TESTS=OFF
cmake --build build -j $(nproc)
sudo cmake --install build          # installs to /usr/local by default
```

If you prefer a user-local install, add `-DCMAKE_INSTALL_PREFIX=$HOME/.local` and set
`CMAKE_PREFIX_PATH` accordingly when building QEPP.

---

## 3 — Install Quantum ESPRESSO (optional, for running calculations)

### From source (recommended — ensures `bands.x` is built)

```bash
git clone https://gitlab.com/QEF/q-e.git qe-src
cd qe-src
cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j $(nproc)
# Add to PATH in ~/.bashrc:
# export QE_HOME=/path/to/qe-src
# export PATH=$QE_HOME/build/bin:$PATH
```

### Via package manager (may lack bands.x)

```bash
sudo apt install quantum-espresso   # Ubuntu — may be outdated
```

---

## 4 — Build QEPP

```bash
git clone https://github.com/albert-hzbn/QEPP.git
cd QEPP
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j $(nproc)
```

The `qepp` binary is placed at `build/qepp`.

### Optional — install to PATH

```bash
sudo cmake --install build          # installs to /usr/local/bin/qepp
```

Or copy manually:

```bash
cp build/qepp ~/.local/bin/         # ensure ~/.local/bin is on PATH
```

---

## 5 — gnuplot at Runtime

Matplot++ requires `gnuplot` to be on `PATH` when generating PNG files.  
If you manage environments with `conda` / `micromamba`, you can install gnuplot there:

```bash
micromamba create -n gnuplot-cf -c conda-forge gnuplot
# Then run qepp through that environment:
micromamba run -n gnuplot-cf qepp band si.bands.dat 6.255 si_bands L,G,X,W,K,G
```

Or export the environment's bin directory into your session PATH:

```bash
export PATH=$(micromamba run -n gnuplot-cf sh -c 'echo $PATH'):$PATH
```

---

## 6 — Quick Verification

```bash
./build/qepp          # prints usage — confirms the build is working
gnuplot --version     # confirms gnuplot is available for plotting
```

Expected output for the first command:

```
Usage:
  qepp cif  <input.cif> <kspacing/Å⁻¹> [output.in] [ecutwfc] [ecutrho]
  qepp dos  <qe.dos> [fermi_eV|qe.out] [prefix]
  qepp band <filband> [fermi_eV|qe.out] [prefix] [L,G,X,...]
  qepp bandqe <input.cif> <scf.in> [bands_pw.in] [bands_pp.in] [pts_per_seg]
```
