# QEPP — Quantum ESPRESSO Pre & Post Processing Tool

A C++17 command-line tool that automates input generation and output analysis for
[Quantum ESPRESSO](https://www.quantum-espresso.org/) DFT calculations.
All plots are rendered as high-resolution PNG files via [Matplot++](https://alandefreitas.github.io/matplotplusplus/).

---

## Features

| Command | Sub-flag | Description |
|---------|----------|-------------|
| `cif`     | `-pre`  | Generate a QE `pw.x` SCF input from a CIF structure file |
| `dos`     | `-post` | Plot total DOS and partial DOS (PDOS) from `dos.x` / `projwfc.x` output |
| `band`    | `-pre`  | Auto-generate band-structure inputs (`pw.x` + `bands.x`) from a CIF |
| `band`    | `-post` | Plot electronic band structure from `bands.x` output |
| `kpath`   | `-pre`  | Suggest and print the standard high-symmetry k-path for a structure |
| `elastic` | `-pre`  | Generate strain-deformed SCF inputs for elastic constant calculation |
| `elastic` | `-post` | Post-process converged strain calculations → full elastic property report |

### cif — SCF Input Generation
- Parses CIF files (fractional coordinates, lattice parameters, symmetry)
- Automatically computes a Monkhorst–Pack k-mesh from a target k-spacing (Å⁻¹)
- Embeds space-group name and IT number as comments
- Configurable cutoff energies (`ecutwfc`, `ecutrho`)

### dos — Density of States Plotting
- Parses QE `.dos` output files (from `dos.x`)
- Shifts energy axis to the Fermi level (accepts a bare number or a QE `.out` file)
- Outputs a `.dos.dat` data file and a `.dos.png` total DOS plot
- **PDOS auto-detection**: if `projwfc.x` partial-DOS files (`prefix.pdos_atm#*`) are present
  in the same directory, two additional plots are generated automatically:
  - `.pdos_elem.png` — elemental contributions (one curve per element + total DOS)
  - `.pdos_orb.png` — orbital-type contributions (s, p, d, f sums + total DOS)
  - Spin-polarised output is detected and handled automatically

### band — Band Structure Plotting
- Parses QE `filband` output (the `nbnd`/`nks` header format produced by `bands.x`)
- Automatically detects high-symmetry segment boundaries from k-vector direction changes
- Draws dashed vertical lines and x-axis tick labels at each high-symmetry point
- Accepts optional comma-separated path labels (e.g. `L,G,X,W,K,G`); `G` renders as Γ
- Fermi energy can be supplied as a number or extracted automatically from a QE output file
- Outputs a `.band.dat` data file and a `.band.png` plot

### kpath — K-Path Suggestion
- Determines the crystal family (cF, cI, cP, …) from the CIF lattice and space group
- Prints the suggested standard k-path and the full `K_POINTS crystal_b` block to stdout
- Optionally writes the block to a file for pasting into a QE input
- Configurable points per segment (default: 20)

### elastic -pre — Strain-Deformed Input Generation
- Detects crystal family automatically via spglib (cubic, hexagonal, tetragonal, …)
- Selects the minimal independent set of strain patterns for the detected symmetry
- Writes one deformed `pw.x` input per (pattern, delta) combination
- Stores only the user-chosen run parameters (`ndeltas`, `max_delta`) in `elastic_setup.dat`;
  all structure-derived quantities (family, volume, patterns) are re-derived at post-processing
  time from the original template — no redundant data to go stale
- Propagates `occupations`, `smearing`, and `degauss` so metallic systems work without
  any manual editing of the deformed inputs

### elastic -post — Elastic Property Report
Post-processes the converged strain calculations and prints a comprehensive report:

**Stiffness and compliance tensors** (full 6×6)

**VRH-averaged elastic moduli** (Voigt/Reuss/Hill):
- Bulk modulus K, shear modulus G, Young's modulus E, Poisson's ratio ν
- Zener anisotropy A_Z

**Lamé constants** λ and μ

**Mechanical character & anisotropy indices:**
- Pugh's ratio K/G (ductile/brittle criterion)
- Cauchy pressure C₁₂ − C₄₄ (bonding character)
- Universal anisotropy index A^U (Ranganathan & Ostoja-Starzewski 2008)
- Percent bulk/shear anisotropy A_B, A_G
- Per-plane shear anisotropy A₁, A₂, A₃ (Chung & Buessem 1967)

**Vickers hardness** — four independent empirical models:
| Model | Formula | Reference |
|-------|---------|-----------|
| Chen 2011 | 2(k²G)^0.585 − 3 | Intermetallics 19 (2011) 1275 |
| Tian 2012 | 0.92 k^1.137 G^0.708 | IJRMHM 33 (2012) 93 |
| Teter 1998 | 0.151 G_V | MRS Bull. 23 (1998) 22 |
| Niu 2012 | G(1−2ν)/3 | J. Phys.: CM 24 (2012) 405401 |

**Directional elastic properties** (Nye 1957):
- E along [100], [010], [001] from diagonal compliance
- E along [110], [111] for cubic systems
- Linear compressibility β_x, β_y, β_z

**Sound velocities & thermal properties:**
- Longitudinal v_L, transverse v_T, mean/Debye v_m
- Debye temperature Θ_D (Anderson 1963)
- Acoustic Grüneisen parameter γ
- Minimum thermal conductivity κ_min (Cahill–Watson–Pohl 1992)
- Lattice thermal conductivity κ_Slack at 300 K (Slack 1973)
- Empirical melting temperature T_m (Fine, Brown & Marcus 1984)

**Born mechanical stability** check (Born & Huang 1954; Mouhat & Coudert 2014)

---

## Usage

```
qepp <command> -pre  [arguments]   (pre-processing)
qepp <command> -post [arguments]   (post-processing)
qepp help [command] [-pre|-post]
```

```
qepp cif     -pre  <input.cif> <kspacing/Å⁻¹> [output.in] [ecutwfc] [ecutrho]
qepp dos     -post <qe.dos> [fermi_eV|qe.out] [prefix]
qepp band    -pre  <input.cif> <scf.in> [bands_pw.in] [bands_pp.in] [pts_per_seg]
qepp band    -post <filband> [fermi_eV|qe.out] [prefix] [L,G,X,...]
qepp kpath   -pre  <input.cif> [pts_per_segment] [output.kpath]
qepp elastic -pre  <scf_template.in> <outdir> [ndeltas] [max_delta]
qepp elastic -post <scf_template.in> <outdir>
```

### Examples

**Generate an SCF input from a CIF:**
```bash
qepp cif -pre silicon.cif 0.15 si_scf.in 40 320
```

**Plot total DOS and PDOS (Fermi energy read from QE output):**
```bash
# Run dos.x and projwfc.x first, then:
qepp dos -post si.dos si_scf.out si_dos_plot
# → si_dos_plot.dos.dat       (data file)
#   si_dos_plot.dos.png       (total DOS)
#   si_dos_plot.pdos_elem.png (elemental PDOS — if projwfc.x files present)
#   si_dos_plot.pdos_orb.png  (orbital PDOS  — if projwfc.x files present)
```

**Plot a band structure with labelled k-path:**
```bash
qepp band -post si.bands.dat si_scf.out si_bands L,G,X,W,K,G
# → si_bands.band.dat  si_bands.band.png
```

**Auto-generate band-path inputs and run the calculation:**
```bash
# 1. Generate inputs
qepp band -pre silicon.cif si_scf.in si_bands.in si_bands_pp.in 30

# 2. Run pw.x (bands)
pw.x < si_bands.in > si_bands.out

# 3. Run bands.x
bands.x < si_bands_pp.in > si_bands_pp.out

# 4. Plot
qepp band -post silicon.bands.dat si_scf.out si_bands L,G,X,W,K,G
```

**Check the k-path for a structure:**
```bash
qepp kpath -pre silicon.cif
# Prints crystal family, space group, path (L-G-X-W-K-G), and K_POINTS crystal_b block

# Save to file with 30 points per segment:
qepp kpath -pre silicon.cif 30 silicon.kpath
```

**Compute elastic constants (example: Si primitive cell):**
```bash
# 1. Generate deformed inputs (5 strain points, ±3%)
qepp elastic -pre si_scf.in results/elastic 5 0.03
# → results/elastic/e1/{p0030,p0015,p0000,m0015,m0030}/si.in
#   results/elastic/e1e2/…  results/elastic/e4/…

# 2. Run all DFT calculations
for d in results/elastic/*/*/; do
  mkdir -p "$d/tmp"
  pw.x < "$d/si.in" > "$d/si.out"
done

# 3. Post-process → full property report
qepp elastic -post si_scf.in results/elastic
# → results/elastic/elastic_results.txt  (printed to stdout as well)
```

---

## Output Files

| File / Extension | Content |
|------------------|---------|
| `.dos.dat` | Two-column (E−E_F, DOS) text data |
| `.dos.png` | Total DOS plot (Matplot++) |
| `.pdos_elem.png` | Elemental PDOS contributions plot |
| `.pdos_orb.png` | Orbital-type PDOS (s/p/d/f) plot |
| `.band.dat` | Two-column (k-dist, E−E_F) text data, bands separated by blank lines |
| `.band.png` | Band structure plot with k-path labels (Matplot++) |
| `elastic_setup.dat` | Run parameters: `ndeltas`, `max_delta` |
| `elastic_results.txt` | Full elastic property report (mirrors stdout) |

---

## Project Structure

```
QEPP/
├── CMakeLists.txt
├── INSTALL.md
├── README.md
├── include/qe/
│   ├── types.hpp       # Shared data structures (ElasticResults, QEInput, …)
│   ├── utils.hpp       # String helpers
│   ├── cif.hpp         # CIF parsing, k-mesh, k-path suggestion
│   ├── qe_input.hpp    # QE input file parsing and writers
│   ├── dos.hpp         # DOS / PDOS parsing and plotting
│   ├── band.hpp        # Band parsing and plotting
│   └── elastic.hpp     # Elastic constants: setup and post-processing
├── src/
│   ├── main.cpp        # CLI entry point (-pre/-post dispatch)
│   ├── utils.cpp
│   ├── cif.cpp
│   ├── qe_input.cpp
│   ├── dos.cpp         # Total DOS + PDOS plotting via Matplot++
│   ├── band.cpp
│   └── elastic.cpp     # Strain pattern generation, energy fitting, VRH, report
├── tests/
│   ├── fixtures/
│   │   ├── cif/        # CIF test structures
│   │   └── qe-input/   # Reference QE input files
│   ├── band_kpath/     # Si band structure reference data
│   ├── elastic_si_prim/
│   │   ├── elastic_setup.dat     # ndeltas=5, max_delta=0.03
│   │   └── elastic_results.txt  # C11=164, C12=66, C44=105 GPa
│   ├── elastic_si_low/
│   └── elastic_al/
└── third_party/        # (not committed — see INSTALL.md)
    └── matplot-install/
```

---

## Dependencies

| Dependency | Notes |
|------------|-------|
| C++17 compiler (GCC ≥9 / Clang ≥10) | |
| CMake ≥ 3.16 | |
| [Eigen3](https://eigen.tuxfamily.org/) | Header-only; `apt install libeigen3-dev` |
| [Matplot++](https://alandefreitas.github.io/matplotplusplus/) | Build from source; see INSTALL.md |
| [gnuplot](http://www.gnuplot.info/) ≥ 5.4 | Runtime renderer for Matplot++; `apt install gnuplot` |
| [spglib](https://spglib.github.io/spglib/) | Optional but recommended; `apt install libspglib-dev` |
| Quantum ESPRESSO ≥ 7.x | `pw.x`, `dos.x`, `bands.x`, `projwfc.x` must be on `PATH` |

See [INSTALL.md](INSTALL.md) for full build instructions.

---

## Validated Reference: Si (diamond cubic)

Input: 2-atom primitive cell, PBE, `ecutwfc=30 Ry`, 4×4×4 k-mesh, ±3 % strain, 5 points.

| Property | Computed | Literature |
|----------|----------|------------|
| C₁₁ | 164.1 GPa | 165–167 GPa |
| C₁₂ | 65.7 GPa | 63–65 GPa |
| C₄₄ | 104.8 GPa | 79–80 GPa (LDA closer) |
| K_H | 98.5 GPa | 97–99 GPa |
| G_H | 77.3 GPa | 66–68 GPa |
| Θ_D | 697 K | 640–680 K |
| Pugh K/G | 1.27 (brittle) | brittle ✓ |
| Cauchy C₁₂−C₄₄ | −39 GPa (covalent) | covalent ✓ |

---

## License


MIT
