# Elastic Strain

Computes the per-atom elastic strain tensor given a crystal structure type and lattice constant.

## Install

```bash
vpm install @voltlabs/elastic-strain
```

## CLI

```bash
elastic-strain <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--clusters_table <path>` | yes | inferred from upstream | Cluster graph table from an upstream structure-identification stage (PTM / ACNA / PSM). |
| `--clusters_transitions <path>` | yes | inferred from upstream | Cluster transitions table from the same upstream stage. |
| `--neighbor_lattice <path>` | yes | inferred from upstream | Per-atom neighbor topology parquet from the same upstream stage. |
| `--lattice_constant <float>` | no | `3.615` | Lattice constant a₀. |
| `--crystal_structure <type>` | no | `FCC` | Crystal structure: `FCC`, `BCC`, or `HCP`. |
| `--ca_ratio <float>` | no | `1.0` | c/a ratio for hexagonal crystals. |
| `--push_forward` | no | `false` | Push strain tensors to the spatial (Euler) frame. |
| `--calc_deformation_gradient` | no | `true` | Compute the deformation gradient F. |
| `--calc_strain_tensors` | no | `true` | Compute strain tensors. |
| `--threads <int>` | no | auto (physical cores) | Max worker threads. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_elastic_strain.parquet` | Elastic Strain | — (listing only) |
| `{output_base}_atoms.parquet` | Elastic Strain Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/elastic-strain
