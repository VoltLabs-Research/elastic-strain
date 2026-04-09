# ElasticStrain

`ElasticStrain` reconstructs elastic strain fields from an upstream crystal-state package.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- ElasticStrain
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
| `--clusters-table <path>` | Yes | Path to `*_clusters.table` exported upstream. | |
| `--clusters-transitions <path>` | Yes | Path to `*_cluster_transitions.table` exported upstream. | |
| `--crystalStructure <type>` | No | Input crystal structure, such as `BCC`, `FCC`, `HCP`, `SC`, `CUBIC_DIAMOND`, `HEX_DIAMOND`. | `BCC` |
| `--latticeConstant <float>` | Yes | Lattice constant `a₀`. | |
| `--caRatio <float>` | No | `c/a` ratio for hexagonal crystals. | `1.0` |
| `--pushForward` | No | Push strain tensors to the spatial frame. | `false` |
| `--calcDeformationGradient` | No | Compute deformation gradient `F`. | `true` |
| `--calcStrainTensors` | No | Compute strain tensors. | `true` |
| `--threads <int>` | No | Maximum worker threads. | auto capped to physical cores |
| `--help` | No | Print CLI help. | |
