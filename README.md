# QEPP — Quantum ESPRESSO Pre & Post Processing Tool

A C++17 command-line tool that automates input generation and output analysis for
[Quantum ESPRESSO](https://www.quantum-espresso.org/) DFT calculations.
All plots are rendered as high-resolution PNG files via [Matplot++](https://alandefreitas.github.io/matplotplusplus/).

---

## Features

| Command | Sub-flag | Description |
|---------|----------|-------------|
| `cif`     | `-pre`  | Generate a QE `pw.x` SCF input from a CIF structure file |
| `dos`     | `-post` | Plot total DOS/PDOS and estimate d-band center/width from `dos.x` / `projwfc.x` output |
| `band`    | `-pre`  | Auto-generate band-structure inputs (`pw.x` + `bands.x`) from a CIF |
| `band`    | `-post` | Plot electronic band structure from `bands.x` output |
| `kpath`   | `-pre`  | Suggest and print the standard high-symmetry k-path for a structure |
| `elastic` | `-pre`  | Generate strain-deformed SCF inputs for elastic constant calculation |
| `elastic` | `-post` | Post-process converged strain calculations → full elastic property report |
| `charge`  | `-pre`  | Generate `pp.x` inputs for charge density, charge difference, and ELF |
| `charge`  | `-post` | Plot volumetric data from Gaussian cube files |
| `mag`     | `-post` | Summarize per-atom magnetic moments from QE `pw.x` output |
| `stm`     | `-pre`  | Generate `pp.x` ILDOS input (`plot_num=5`) for STM simulation |
| `stm`     | `-post` | Plot 2D constant-height STM maps from cube files |
| `bader`   | `-post` | Parse Henkelman `ACF.dat` and report per-atom Bader charges |
| `conv`    | `-pre`  | Generate convergence sweep SCF inputs for `ecutwfc` or `kspacing` |
| `conv`    | `-post` | Parse sweep energies and produce convergence table + plot |
| `struct`  | `-post` | Print structure summary and estimate Warren-Cowley SRO from QE input/output |
| `parse`   | `-post` | Parse QE `pw.x` output summary (energy, Fermi, force, pressure, timing) |
| `phonon`  | `-pre`  | Generate all DFPT phonon inputs (`ph.x`, `q2r.x`, `matdyn.x`) from an SCF input |
| `phonon`  | `-post` | Post-process matdyn.x outputs → DOS plot, band plot, harmonic thermodynamics |
| `phonon`  | `-dos`  | Plot phonon DOS (phonopy `total_dos.dat` or matdyn.x format) |
| `phonon`  | `-band` | Plot phonon band structure (phonopy `band.yaml` or matdyn.x `.freq` format) |
| `phonon`  | `-ha`   | Harmonic approximation thermodynamics: ZPE, F_vib(T), S(T), Cv(T) |
| `qha`     | `-pre`  | Generate N volume-scaled SCF inputs for a quasi-harmonic approximation grid |
| `qha`     | `-post` | Fit Birch-Murnaghan EOS and compute V(T), α(T), B_T(T), Cp(T), γ(T) |

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
- **d-band descriptors** (when d-PDOS exists): writes `.dband.txt` with
  - **old methodology**: occupied-only moments ($E \le E_F$)
  - **new methodology**: full d-band moments (all sampled energies)
  - reported values: d-band center and d-band width ($\sigma$ from second central moment)

  Using the d-projected DOS $D_d(E)$, QEPP reports:
  - $\mu_d = \dfrac{\int E\,D_d(E)\,dE}{\int D_d(E)\,dE}$
  - $\sigma_d = \sqrt{\dfrac{\int E^2\,D_d(E)\,dE}{\int D_d(E)\,dE} - \mu_d^2}$

  where the integration range depends on the selected methodology:
  - old: up to $E_F$
  - new: full sampled energy range

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
Post-processes the converged strain calculations and computes the full elastic tensor,
VRH-averaged moduli (K, G, E, ν), Lamé constants, anisotropy indices, Vickers hardness
(four empirical models), sound velocities, Debye temperature, and Born mechanical stability check.

### conv — Convergence Sweep and Analysis
- `conv -pre` creates a parameter sweep over either `ecutwfc` (Ry) or `kspacing` (1/Ang)
  from a template SCF input and writes one case per subdirectory.
- `conv -post` scans completed runs, parses total energies, computes
  $|\Delta E|$ in meV/atom relative to the most converged point, and writes:
  - `<prefix>.conv.txt` (table)
  - `<prefix>.conv.png` (convergence plot)

### struct — Structure Summary from QE Input
- Parses a QE SCF input and prints lattice constants, angles, cell volume,
  lattice vectors, and per-atom fractional/Cartesian coordinates.
- Supports explicit `CELL_PARAMETERS` (`ibrav=0`) and common cubic
  `ibrav` fallbacks (`1`, `2`, `3`) using `celldm(1)`/`A`.
- Writes `<prefix>.struct.txt`.
- Optional Warren-Cowley SRO analysis (`--sro`) for neighbor shells:
  - $\alpha_{ij}^{(s)} = 1 - P(j|i,s)/c_j$
  - supports both QE input and QE output structures (`--source input|output|auto`)
  - writes `<prefix>.sro.txt`

### phonon — DFPT Phonon Workflow

A complete interface to Quantum ESPRESSO's DFPT machinery via `ph.x`, `q2r.x`, and `matdyn.x`.

#### phonon -pre — Input Generation
- Reads an existing pw.x SCF input and writes four ready-to-run QE input files:
  - `ph.in` — `ph.x` phonon calculation (ldisp mode, q-mesh)
  - `q2r.in` — `q2r.x` real-space interatomic force constant conversion
  - `matdyn_dos.in` — `matdyn.x` phonon density of states (finer q-mesh)
  - `matdyn_band.in` — `matdyn.x` phonon dispersion along high-symmetry path
- Configurable q-mesh (`--nq`), DOS mesh (`--nq-dos`), convergence threshold, acoustic sum rule, and Born effective charges flag

#### phonon -post — Unified Post-Processing
Convenience wrapper that runs DOS, band-structure, and/or HA analysis in one call from matdyn.x outputs:
```
qepp phonon -post <prefix> [--dos] [--band] [--ha] [--natom N] [--labels G,X,...] [--tmax T]
```
If no flags are given, all three analyses are performed.

#### phonon -dos / -band — Individual Plots
- **-dos**: Parses matdyn.x or phonopy `total_dos.dat`; plots total and (optionally) projected DOS with configurable atom labels
- **-band**: Parses matdyn.x `.freq` or phonopy `band.yaml`; draws dispersion with labelled high-symmetry points

#### phonon -ha — Harmonic Approximation Thermodynamics
Integrates the phonon DOS to yield:
$$F_{\rm vib}(T) = k_BT\int_0^\infty g(\nu)\ln\!\left(2\sinh\tfrac{h\nu}{2k_BT}\right)d\nu$$
$$C_V(T) = k_B\int_0^\infty g(\nu)\left(\frac{h\nu}{k_BT}\right)^2\!\frac{e^{h\nu/k_BT}}{(e^{h\nu/k_BT}-1)^2}\,d\nu$$
- Reports ZPE, F_vib(T), S(T), Cv(T), E_vib(T) per atom
- Both matdyn.x (cm⁻¹) and phonopy (THz) formats are auto-detected

---

### qha — Quasi-Harmonic Approximation

Computational workflow to obtain finite-temperature thermodynamic properties by
combining electronic energies with phonon free energies across a volume grid.

#### qha -pre — Volume Grid Generation
- Reads a pw.x SCF input and generates N symmetrically scaled copies
  spanning ±range% of the equilibrium volume
- Writes `<outDir>/v01/ … <outDir>/vNN/` each with a ready-to-run SCF input
- Writes a `qha_summary.in` template to be filled in after running SCF + DFPT at each volume

#### qha -post — Thermal Property Calculation
Reads `qha_summary.in` (columns: volume Å³, energy Ry, path to matdyn DOS), then:
1. Computes HA vibrational free energy at each volume over the full T range
2. Fits a Birch-Murnaghan 3rd-order EOS to E_static(V) and to G(V,T) at each T
3. Derives thermal properties from the free-energy minimum:
   - V(T), V/atom(T)
   - Isothermal bulk modulus $B_T(T)$
   - Volumetric thermal expansion $\alpha(T) = \frac{1}{V}\frac{dV}{dT}$
   - Isochoric and isobaric heat capacities $C_V$, $C_P$
   - Vibrational entropy S(T)
   - Macroscopic Grüneisen parameter $\gamma = V\alpha B_T/C_V$

---

### phonon — DFPT Phonon Workflow

A complete interface to Quantum ESPRESSO's DFPT machinery via `ph.x`, `q2r.x`, and `matdyn.x`.

#### phonon -pre — Input Generation
- Reads an existing pw.x SCF input and writes four ready-to-run QE input files:
  - `ph.in` — `ph.x` phonon calculation (ldisp mode, q-mesh)
  - `q2r.in` — `q2r.x` real-space interatomic force constant conversion
  - `matdyn_dos.in` — `matdyn.x` phonon density of states (finer q-mesh)
  - `matdyn_band.in` — `matdyn.x` phonon dispersion along high-symmetry path
- Configurable q-mesh (`--nq`), DOS mesh (`--nq-dos`), convergence threshold, acoustic sum rule, and Born effective charges flag (`--epsil`)

#### phonon -post — Unified Post-Processing
Convenience wrapper that runs DOS, band-structure, and/or HA analysis from matdyn.x outputs:
- Reads `<prefix>.phonon.dos` and/or `<prefix>.phonon_band.freq`
- If no `--dos/--band/--ha` flag is given, all three analyses are performed

#### phonon -dos / -band — Individual Plots
- **-dos**: Parses matdyn.x or phonopy `total_dos.dat`; plots total DOS (and optionally projected DOS with configurable atom labels)
- **-band**: Parses matdyn.x `.freq` or phonopy `band.yaml`; draws dispersion with labelled high-symmetry points (`G` renders as Γ)

#### phonon -ha — Harmonic Approximation Thermodynamics
Integrates the phonon DOS to yield vibrational free energy, entropy, and heat capacity:

$$C_V(T) = k_B\int_0^\infty g(\nu)\!\left(\frac{h\nu}{k_BT}\right)^{\!2}\frac{e^{h\nu/k_BT}}{(e^{h\nu/k_BT}-1)^2}\,d\nu$$

- Reports ZPE, F_vib(T), S(T), Cv(T), E_vib(T) per atom and per unit cell
- Both matdyn.x (cm⁻¹, auto-detected) and phonopy (THz) formats are supported

---

### qha — Quasi-Harmonic Approximation

Computes finite-temperature thermodynamic properties by combining static DFT energies
with phonon free energies across a volume grid.

#### qha -pre — Volume Grid Generation
- Reads a pw.x SCF input and generates N symmetrically scaled copies spanning ±range% of the equilibrium volume
- Writes `<outDir>/v01/ … <outDir>/vNN/` each with a ready-to-run SCF input
- Writes a `qha_summary.in` template (volume, energy placeholder, DOS path) to be completed after running SCF + DFPT at each volume

#### qha -post — Thermal Property Calculation
Reads `qha_summary.in` (columns: volume Å³, energy Ry, path to matdyn DOS), then:
1. Computes HA vibrational free energy at each volume over the full T range
2. Fits a 3rd-order Birch-Murnaghan EOS to E_static(V) and to G(V,T) at each T
3. Derives thermal properties from the free-energy minimum at each T:

| Output | Symbol |
|--------|--------|
| Equilibrium volume | V(T) |
| Isothermal bulk modulus | $B_T(T)$ |
| Volumetric thermal expansion | $\alpha(T) = \frac{1}{V}\frac{dV}{dT}$ |
| Isochoric / isobaric heat capacity | $C_V$, $C_P = C_V + TV\alpha^2 B_T$ |
| Vibrational entropy | S(T) |
| Macroscopic Grüneisen parameter | $\gamma = V\alpha B_T / C_V$ |

---

### parse — QE Output Summary Parser
- Parses a QE `pw.x` output file and extracts:
  - Total energy (Ry/eV)
  - Fermi level
  - Per-atom forces and max force
  - Pressure
  - SCF convergence + iteration count
  - Wall time
- Writes `<prefix>.parse.txt`.

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
qepp charge  -pre  <scf.in> [outdir]
qepp charge  -post <cube_file> [prefix]
qepp mag     -post <qe.out> [prefix]
qepp stm     -pre  <scf.in> [bias_eV] [outdir]
qepp stm     -post <stm.cube> [prefix] [height_ang]
qepp bader   -post <ACF.dat> [scf.in] [prefix]
qepp conv    -pre  <scf.in> <ecutwfc|kspacing> <min> <max> <step> [outdir]
qepp conv    -post <outdir> <ecutwfc|kspacing> [prefix]
qepp struct  -post <structure_file> [prefix] [--sro] [--nshells N] [--tol T] [--source input|output|auto]
qepp parse   -post <qe.out> [prefix]
qepp phonon  -pre  <scf.in> [outdir] [--nq NQ1 NQ2 NQ3] [--nq-dos NQ1 NQ2 NQ3] [--epsil]
qepp phonon  -post <prefix> [outprefix] [--dos] [--band] [--ha] [--natom N] [--labels L,...] [--tmin T] [--tmax T] [--dt T]
qepp phonon  -dos  <dos_file> [prefix] [--pdos projected_dos.dat] [--labels L,...]
qepp phonon  -band <freq_or_yaml> [prefix] [--labels G,X,W,L,G]
qepp phonon  -ha   <dos_file> [prefix] [--natom N] [--tmin T] [--tmax T] [--dt T]
qepp qha     -pre  <scf.in> [--nvolumes N] [--range R] [--outdir D]
qepp qha     -post <qha_summary.in> [output_prefix] [--tmin T] [--tmax T] [--dt T]
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
#   si_dos_plot.dband.txt     (d-band center/width report — if d-PDOS exists)
```

**Inspect d-band descriptor output:**
```text
=== d-band descriptors from PDOS ===
Method: old (occupied states up to E_F)
  d-band center (eV)    : -1.234567
  center - E_F (eV)     : -2.345678
  d-band width (eV)     : 1.456789
  integrated d-DOS area : 5.432100

Method: new (full d-band moment over all states)
  d-band center (eV)    : -0.987654
  center - E_F (eV)     : -2.098765
  d-band width (eV)     : 1.678901
  integrated d-DOS area : 8.765400
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

**Generate and analyze an `ecutwfc` convergence sweep:**
```bash
# 1. Generate sweep inputs
qepp conv -pre si_scf.in ecutwfc 20 80 10 results/conv_ecut

# 2. Run QE in each case
for d in results/conv_ecut/*/; do
  pw.x < "$d/scf.in" > "$d/scf.out"
done

# 3. Summarize convergence
qepp conv -post results/conv_ecut ecutwfc results/conv_ecut/ecutwfc
# → results/conv_ecut/ecutwfc.conv.txt
#   results/conv_ecut/ecutwfc.conv.png
```

**Print structure info from an SCF input:**
```bash
qepp struct -post si_scf.in si_struct
# → si_struct.struct.txt
```

**Estimate Warren-Cowley SRO from input/output structure:**
```bash
# From QE input
qepp struct -post ni_co.scf.in ni_co --sro --nshells 3
# From QE relax/vc-relax output (last complete final block)
qepp struct -post ni_co.relax.out ni_co --sro --source output
# → ni_co.struct.txt
#   ni_co.sro.txt
```

**Parse QE SCF output summary:**
```bash
qepp parse -post si_scf.out si_summary
# → si_summary.parse.txt
```

**Run the full DFPT phonon workflow for Si:**
```bash
# 1. Generate all DFPT inputs from SCF input
qepp phonon -pre si_scf.in phonon_inputs/ --nq 4 4 4 --nq-dos 16 16 16

# 2. Run QE calculations
cd phonon_inputs/
mpirun -np 4 pw.x  -input si_scf.in > si_scf.out
mpirun -np 4 ph.x  -input ph.in     > ph.out
q2r.x               < q2r.in          > q2r.out
matdyn.x            < matdyn_dos.in   > matdyn_dos.out
matdyn.x            < matdyn_band.in  > matdyn_band.out

# 3. Post-process all at once
qepp phonon -post si --natom 2 --labels G,X,W,L,G --tmax 1500
# → si.phonon_dos.png   si.phonon_band.png   si.phonon_ha.txt   si.phonon_ha.png
```

**Compute QHA thermal properties (quasi-harmonic approximation):**
```bash
# 1. Generate 9-volume grid spanning ±6%
qepp qha -pre si_scf.in --nvolumes 9 --range 12 --outdir si_qha
# → si_qha/v01/ … si_qha/v09/  (scaled SCF inputs)
#   si_qha/qha_summary.in       (template; fill in energies after SCF runs)

# 2. Run SCF + DFPT at each volume (parallel, e.g. with a loop)
for d in si_qha/v*/; do
  cd "$d" && mpirun -np 4 pw.x -input si_scf.in > qe.out
  # run ph.x → q2r.x → matdyn.x (writes si.phonon.dos)
  cd -
done

# 3. Fill energies in qha_summary.in, then post-process
qepp qha -post si_qha/qha_summary.in si_qha_results --tmin 0 --tmax 1500 --dt 10
# → si_qha_results.qha.txt  (V(T), α(T), B_T(T), Cv(T), Cp(T), S(T), γ(T))
```

---

## Output Files

| File / Extension | Content |
|------------------|---------|
| `.dos.dat` | Two-column (E−E_F, DOS) text data |
| `.dos.png` | Total DOS plot (Matplot++) |
| `.pdos_elem.png` | Elemental PDOS contributions plot |
| `.pdos_orb.png` | Orbital-type PDOS (s/p/d/f) plot |
| `.dband.txt` | d-band center and width (old occupied-only + new full-band methods) |
| `.band.dat` | Two-column (k-dist, E−E_F) text data, bands separated by blank lines |
| `.band.png` | Band structure plot with k-path labels (Matplot++) |
| `elastic_setup.dat` | Run parameters: `ndeltas`, `max_delta` |
| `elastic_results.txt` | Full elastic property report (mirrors stdout) |
| `.mag.txt` | Per-atom magnetic moments + total/absolute magnetization summary |
| `.stm.png` | Constant-height STM 2D map |
| `.bader.txt` | Per-atom Bader charge/volume summary |
| `.conv.txt` | Convergence table for sweep results |
| `.conv.png` | Convergence curve plot ($|\Delta E|$ in meV/atom) |
| `.struct.txt` | Structure summary from QE SCF input |
| `.sro.txt` | Warren-Cowley SRO table by neighbor shell and element pair |
| `.parse.txt` | Parsed QE output summary |
| `.phonon_dos.dat` | Two-column (freq THz, DOS) phonon density of states data |
| `.phonon_dos.png` | Phonon DOS plot |
| `.phonon_band.dat` | Phonon band structure data (k-distance, freq THz per band) |
| `.phonon_band.png` | Phonon dispersion plot with high-symmetry labels |
| `.phonon_ha.txt` | HA thermodynamics table: T, F_vib, S, Cv, E_vib (per atom) |
| `.phonon_ha.png` | Three-panel plot: Cv(T), S(T), F_vib(T) |
| `.qha.txt` | QHA thermal properties: T, V, α, B_T, Cv, Cp, S, γ |

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

## License


MIT
