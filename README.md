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
| `qha`         | `-pre`  | Generate N volume-scaled SCF inputs for a quasi-harmonic approximation grid |
| `qha`         | `-post` | Fit Birch-Murnaghan EOS and compute V(T), α(T), B_T(T), Cp(T), γ(T) |
| `qha_elastic` | `-pre`  | Generate volume-grid SCF + elastic strain + DFPT inputs for T-dependent elastic constants |
| `qha_elastic` | `-post` | Compute C_ij(T) from quasi-static + phonon bulk-modulus contributions; VRH averages |

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

### elastic -run — Built-In 0 K Workflow Execution
- Runs the equilibrium SCF template and every strained SCF generated by `elastic -pre`
- Skips stages that already contain `JOB DONE`, so reruns resume cleanly
- Supports direct QE launch tuning through `--np`, `--ni`, `--nk`, `--nb`, `--nt`, and `--nd`

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

### qha_elastic — Temperature-Dependent Elastic Constants

Extends the QHA approach to elastic properties by combining the energy-strain method
with vibrational free energies across a volume grid.  Two physical contributions are included:

1. **Quasi-static** (dominant): $C_{ij}^{QS}(T)$ — the static elastic tensor cubic-polynomial-interpolated to the QHA equilibrium volume $V(T)$.
2. **Phonon bulk-modulus correction**: $B^{\rm ph}(T) = V(T)\cdot\partial^2 F_{\rm vib}/\partial V^2\big|_{V(T),T}$ — obtained by fitting a cubic polynomial to the phonon free energy over the volume grid and evaluating the second derivative at $V(T)$.  Added to all normal–normal Voigt components ($i,j \in \{0,1,2\}$).

The total tensor $C_{ij}(T) = C_{ij}^{QS}(T) + \Delta C_{ij}^{\rm ph}(T)$ is then used to compute
VRH-averaged bulk modulus $K_H$, shear modulus $G_H$, Young's modulus $E_H$, and Poisson's ratio $\nu_H$.

#### qha_elastic -pre — Input Generation
- Reads a pw.x SCF input, creates N volume-scaled directories (`v01/ … vNN/`).
- In each volume directory:
  - Writes an isotropically scaled pw.x SCF input.
  - Calls the same strain-pattern engine as `elastic -pre` to write deformed inputs in `elastic/`.
  - Writes DFPT phonon inputs in `dfpt/` (identical to `phonon -pre`).
- Writes `qha_elastic_summary.in` with columns: `volume(Å³)  energy(Ry)  scf_template  elastic_dir  phonon_dos_path`.

Options: `--nvolumes` (default 7), `--range` (percent, default 10), `--outdir`,
`--ndeltas` (strain points per pattern, default 5), `--maxdelta` (max strain, default 0.03).

#### qha_elastic -run — Built-In Workflow Execution
- Runs the top-level SCF, all strain SCFs, and the DFPT chain (`ph.x`, `q2r.x`, `matdyn.x`) for every `v*/` directory.
- Skips stages that already contain `JOB DONE`, so interrupted datasets can be resumed with the same command.
- Accepts QE launch controls `--np`, `--ni`, `--nk`, `--nb`, `--nt`, and `--nd`, plus `--exclude v04,v07` to skip known bad volumes.

#### qha_elastic -post — C_ij(T) Calculation
Reads `qha_elastic_summary.in`, then:
1. Calls `elastic -post` logic for each volume → $C_{ij}^{static}(V_n)$.
2. Calls `phonon -ha` DOS integration for each volume at each T → $F_{\rm vib}(V_n,T)$.
3. Fits cubic polynomials in $V$ to each of the 36 $C_{ij}$ elements and to $F_{\rm vib}$.
4. Runs the full QHA minimisation → $V(T)$, $\alpha(T)$, $B_T^{QS}(T)$.
5. Evaluates $C_{ij}^{QS}(V(T))$ and $B^{\rm ph}(T) = V(T)\cdot d^2F_{\rm vib}/dV^2$ at $V(T)$.
6. Adds the phonon correction to normal–normal Voigt entries and computes VRH averages.

### Reproducing The Al qha_elastic litfix_new Run

The aluminum `qha_elastic` workflow can now be executed entirely through built-in `qepp qha_elastic` commands. No external shell helper is required.

#### 1. Prepare An Equilibrium SCF Template

Create a cubic Al SCF template, for example from the near-minimum volume used here:

```bash
cat > /tmp/al_eq_scf.in <<'EOF'
&CONTROL
  calculation = 'scf'
  prefix      = 'allit'
  outdir      = './tmp'
  pseudo_dir  = '/home/albert/Softwares/qe-7.5/pseudo'
  tprnfor     = .false.
  tstress     = .false.
/
&SYSTEM
  ibrav       = 0,
  nat         = 1,
  ntyp        = 1,
  ecutwfc     = 55.0
  ecutrho     = 440.0
  nosym       = .true.
  noinv       = .true.
  occupations = 'smearing'
  smearing    = 'mp'
  degauss     = 0.02
/
&ELECTRONS
  conv_thr    = 1.0d-10
  mixing_beta = 0.30
/
ATOMIC_SPECIES
  Al 26.9815385 Al.pbe-nl-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Al 0.0 0.0 0.0
K_POINTS automatic
  12 12 12 0 0 0
CELL_PARAMETERS angstrom
        0.0000000000  2.0195947927  2.0195947927
        2.0195947927  0.0000000000  2.0195947927
        2.0195947927  2.0195947927  0.0000000000
EOF
```

#### 2. Generate The 9-Volume Dataset

```bash
qepp qha_elastic -pre /tmp/al_eq_scf.in \
  --nvolumes 9 \
  --range 12 \
  --ndeltas 7 \
  --maxdelta 0.015 \
  --nq 2 \
  --nq-dos 16 \
  --outdir tests/al_qha_el_litfix_new
```

This creates `v01` through `v09`, each with one SCF input, `elastic/` strain inputs, `dfpt/` phonon inputs, and a `qha_elastic_summary.in` file. The `--nq 2` flag writes a `2x2x2` DFPT mesh directly into every generated `ph.in` file, so no manual editing is needed.

#### 3. Run The Full Calculation

```bash
qepp qha_elastic -run tests/al_qha_el_litfix_new --np 20 --nk 4 --nd 5
```

What `qha_elastic -run` does automatically:
- runs the top-level SCF for each volume and writes `qe.out`
- runs all 21 elastic strain calculations in `elastic/*/*/`
- runs `ph.x`, `q2r.x`, and `matdyn.x`
- creates `tmp/_ph0` before `ph.x`
- resumes cleanly by skipping stages that already contain `JOB DONE`
- forwards `--ni`, `--nk`, `--nb`, `--nt`, and `--nd` to the QE executables so the run can be tuned for the target machine

#### 4. Drop v04 And Post-Process

In this run, `v04` failed in `ph.x` with `FFT grid incompatible with symmetry`, while the other 8 volumes completed. The built-in `-post` mode can exclude failed volumes directly and auto-read the static SCF energies from `qe.out`, so the summary file does not need manual editing:

```bash
qepp qha_elastic -post tests/al_qha_el_litfix_new --exclude v04 --tmin 0 --tmax 1000 --dt 100
```

This writes `qha_elastic_summary.qha_elastic.txt` in the dataset directory.

#### 5. Reduced-Set 300 K Result Vs Ab Initio Reference

Using the 8 completed volumes (`v01`, `v02`, `v03`, `v05`, `v06`, `v07`, `v08`, `v09`), the reduced-set `qha_elastic` result at 300 K was:

| Quantity | QEPP reduced set | Ab initio reference | Delta | Relative error |
|---------|------------------|---------------------|-------|----------------|
| C11 | 89.7791 GPa | 107.0 GPa | -17.2209 GPa | -16.09% |
| C12 | 61.9274 GPa | 60.0 GPa | +1.9274 GPa | +3.21% |
| C44 | 35.6840 GPa | 32.0 GPa | +3.6840 GPa | +11.51% |

These values came from `tests/al_qha_el_litfix_new/qha_elastic_summary.qha_elastic.txt` after running `qepp qha_elastic -post tests/al_qha_el_litfix_new --exclude v04 ...`.

Outputs a full column-data report and a concise summary table to stdout.

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
qepp <command> -run  [arguments]   (workflow execution)
qepp <command> -post [arguments]   (post-processing)
qepp help [command] [-pre|-run|-post]
```

```
qepp cif     -pre  <input.cif> <kspacing/Å⁻¹> [output.in] [ecutwfc] [ecutrho]
qepp dos     -post <qe.dos> [fermi_eV|qe.out] [prefix]
qepp band    -pre  <input.cif> <scf.in> [bands_pw.in] [bands_pp.in] [pts_per_seg]
qepp band    -post <filband> [fermi_eV|qe.out] [prefix] [L,G,X,...]
qepp kpath   -pre  <input.cif> [pts_per_segment] [output.kpath]
qepp elastic -pre  <scf_template.in> [outdir] [ndeltas] [max_delta] [--outdir D] [--ndeltas N] [--maxdelta D]
qepp elastic -run  <scf_template.in> <outdir> [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N]
qepp elastic -run  <outdir_or_volume_dir> [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N]
qepp elastic -post <scf_template.in> <outdir>
qepp elastic -post <outdir_or_volume_dir>
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
qepp qha_elastic -pre  <scf.in> [outdir] [--nvolumes N] [--range R] [--outdir D] [--ndeltas N] [--maxdelta D]
qepp qha_elastic -run  <dataset_dir> [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N] [--exclude v04,v07]
qepp qha_elastic -run  [--np N] [--ni N] [--nk N] [--nb N] [--nt N] [--nd N] [--exclude v04,v07]
qepp qha_elastic -post <qha_elastic_summary.in|dataset_dir> [output_prefix] [--tmin T] [--tmax T] [--dt T] [--exclude v04]
qepp qha_elastic -post [--tmin T] [--tmax T] [--dt T] [--exclude v04]
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
# (or rely on the default folder: si_scf_elastic)
qepp elastic -pre si_scf.in
# → results/elastic/e1/{p0030,p0015,p0000,m0015,m0030}/si.in
#   results/elastic/e1e2/…  results/elastic/e4/…

# 2. Run the equilibrium and strained SCFs
qepp elastic -run si_scf.in results/elastic --np 16 --nk 4
# (or if folder names are obvious)
qepp elastic -run results/elastic --np 16 --nk 4

# 3. Post-process → full property report
qepp elastic -post si_scf.in results/elastic
qepp elastic -post results/elastic
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

# 2. Run SCF + DFPT at each volume
# All four executables must run from the volume directory (paths in the
# generated inputs are relative to that directory, not to dfpt/).
for d in si_qha/v*/; do
  cd "$d"
  mpirun -np 4 pw.x  -input si_scf.in          > qe.out              # SCF
  mpirun -np 4 ph.x  -input dfpt/ph.in          > dfpt/ph.out         # DFPT phonons
  q2r.x              < dfpt/q2r.in               > dfpt/q2r.out        # IFCs
  matdyn.x           < dfpt/matdyn_dos.in        > dfpt/matdyn_dos.out # phonon DOS
  # → writes <prefix>.phonon.dos in this directory
  cd -
done

# 2.5 Fill DFT total energies (Ry) into qha_summary.in
# QE prints the converged energy on lines starting with '!'; sed replaces
# each TODO_fill_energy placeholder in document order (one per iteration).
for d in si_qha/v*/; do
  e=$(grep "!.*total energy" "$d/qe.out" | tail -1 | awk '{print $5}')
  sed -i "0,/TODO_fill_energy/s/TODO_fill_energy/${e}/" si_qha/qha_summary.in
done

# 3. Post-process
qepp qha -post si_qha/qha_summary.in si_qha_results --tmin 0 --tmax 1500 --dt 10
# → si_qha_results.qha.txt  (V(T), α(T), B_T(T), Cv(T), Cp(T), S(T), γ(T))
```

**Compute temperature-dependent elastic constants (QHA elastic):**
```bash
# 1. Generate volume grid with elastic strain + DFPT inputs
qepp qha_elastic -pre si_scf.in --nvolumes 7 --range 10 --outdir si_qha_el
# (equivalent positional outdir form)
qepp qha_elastic -pre si_scf.in si_qha_el --nvolumes 7 --range 10
# → si_qha_el/v01/ … si_qha_el/v07/  (each: si.in, elastic/, dfpt/)
#   si_qha_el/qha_elastic_summary.in  (fill energies after SCF runs)

# 2. Run the full workflow with explicit QE parallelization tags
qepp qha_elastic -run si_qha_el --np 20 --nk 4 --nd 5
# (or from inside si_qha_el/, just run)
qepp qha_elastic -run --np 20 --nk 4 --nd 5

# 3. Post-process → C_ij(T) report
qepp qha_elastic -post si_qha_el si_qha_el_results --tmin 0 --tmax 1200 --dt 10
# (or from inside si_qha_el/, using the auto-detected summary)
qepp qha_elastic -post --tmin 0 --tmax 1200 --dt 10
# → si_qha_el_results.qha_elastic.txt
#   Columns: T(K), V(Å³), α(1/K), C11…C66(GPa) [21 upper-triangle],
#            KH, GH, EH, νH (GPa/GPa/GPa/–), B_QS(GPa), B_ph(GPa)
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
| `.qha_elastic.txt` | T-dependent elastic constants: T, V, α, C11–C66 (21 upper-triangle), KH, GH, EH, νH, B_QS, B_ph |

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
