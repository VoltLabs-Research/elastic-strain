# ElasticStrain

`ElasticStrain` reconstructs elastic strain fields from an upstream crystal-state package.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- ElasticStrain
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)
- `structure-identification/1.0.0` (from the `StructureIdentification` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Run

```bash
./build/build/Release/elastic-strain --help
```

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

## CLI

Usage:

```bash
elastic-strain <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--clusters_table <path>` | Yes | Path to `*_clusters.table` exported upstream. | |
| `--clusters_transitions <path>` | Yes | Path to `*_cluster_transitions.table` exported upstream. | |
| `--crystal_structure <type>` | No | Input crystal structure, such as `BCC`, `FCC`, `HCP`, `SC`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. | `BCC` |
| `--lattice_constant <float>` | Yes | Lattice constant `a₀`. | |
| `--ca_ratio <float>` | No | `c/a` ratio for hexagonal crystals. | `1.0` |
| `--push_forward` | No | Push strain tensors to the spatial frame. | `false` |
| `--calc_deformation_gradient` | No | Compute deformation gradient `F`. | `true` |
| `--calc_strain_tensors` | No | Compute strain tensors. | `true` |
| `--threads <int>` | No | Maximum worker threads. | auto capped to physical cores |
| `--help` | No | Print CLI help. | |
