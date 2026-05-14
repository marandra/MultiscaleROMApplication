# MultiscaleROMApplication

A custom application of the [Kratos
Multiphysics](https://github.com/KratosMultiphysics/Kratos) FEM
framework implementing **HPR-FE²** (High-Performance Reduced Finite
Element²), a hyper-reduced order modelling technique for multiscale
simulation of materials with heterogeneous microstructure.

## What it does

The application replaces the inner microscopic FE problem of a
classical FE² scheme by a hyper-reduced surrogate, combining:

- **ROM** — the RVE displacement field is projected onto a
  low-dimensional basis obtained from an SVD of strain snapshots
  collected during an *offline* training phase.
- **Hyper-reduction** — the internal force and tangent stiffness
  integrals over the RVE are evaluated only at a small set of
  optimally selected sampling points, instead of a full loop over the
  Gauss points.

This yields up to **4 orders of magnitude speed-up** over a standard
FE² solver while keeping the accuracy loss below **1%**, making
two-scale simulations of industrial-scale problems tractable.

## Place in the workflow

The application is the solver-side component of a two-piece system:

```
┌──────────────────────────────┐      ┌──────────────────────────────┐
│  hprfe2 (Python)             │      │  MultiscaleROMApplication    │
│                              │      │  (C++/Python, this repo)     │
│  - case generation           │      │                              │
│  - offline SVD training      │ ───► │  - reads trained ROM         │
│  - Slurm job orchestration   │      │  - integrates RVE response   │
│  - result collection         │      │    at sampling points        │
│  - field reconstruction      │ ◄─── │  - couples to macro FE       │
└──────────────────────────────┘      └──────────────────────────────┘
```

`hprfe2` (the orchestration package) lives in
[its own repository](https://github.com/marandra/hprfe2).

## Repository layout

```
applications/MultiscaleROMApplication/
├── custom_processes/        # offline snapshot acquisition,
│                            # SVD-basis import, hyper-reduction
│                            # sampling-point selection
├── custom_strategies/       # online ROM solution strategies
├── custom_elements/         # reduced RVE elements
├── custom_io/               # HDF5 import/export of bases and snapshots
├── python_scripts/          # high-level Python entry points
│                            # (MainKratos.py, OfflineKratos.py, ...)
└── tests/                   # regression tests
```

## Usage

### Offline (training) stage

Run a set of microscopic simulations to collect strain/stress snapshots,
then build the reduced basis and select the hyper-reduction sampling
points:

```bash
python3 python_scripts/MainKratos.py            # snapshot acquisition
mkdir offline_data
python3 python_scripts/OfflineKratos_1.py ..    # SVD + sampling
```

The trained ROM (basis + sampling weights) is written to HDF5 for the
online stage to consume.

### Online (macroscopic) stage

A two-scale macroscopic analysis runs as a standard Kratos simulation,
with `MultiscaleROMApplication` providing the reduced RVE response at
each macroscopic Gauss point. See `python_scripts/` for the entry
points and example problem types.

### Field reconstruction

Full-field RVE responses at user-selected macroscopic points can be
reconstructed from the reduced solution in post-processing; this is
typically handled by the companion `hprfe2` package.

## Build

`MultiscaleROMApplication` follows the standard Kratos application
build flow: add the application to the Kratos applications list and
compile against a Kratos source tree.

```bash
# In your Kratos applications config:
add_app /path/to/MultiscaleROMApplication
```

Refer to the
[Kratos build documentation](https://github.com/KratosMultiphysics/Kratos/wiki)
for the surrounding setup.

## Citing

If you use this code in academic work, please cite:

> M. Raschi, O. Lloberas, A. Huespe, J. Oliver.
> *High performance technique for multiscale finite element (HPR-FE²):
> towards industrial multiscale FE software.*
> Computer Methods in Applied Mechanics and Engineering, 2021.

> M. Caicedo, J. L. Mroginski, S. Toro, M. Raschi, A. Huespe, J. Oliver.
> *High performance reduced order modeling techniques based on optimal
> energy quadrature: application to geometrically non-linear multiscale
> inelastic material modeling.*
> Archives of Computational Methods in Engineering, 2019.

## License

Inherits the Kratos Multiphysics licence (BSD).
