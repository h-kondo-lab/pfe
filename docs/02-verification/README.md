---
navigation:
  order: 200
---

# Verification

Every numerical claim in the paper's tables can be reproduced from the two scripts in [`paper1/`](../../paper1/). No experimental value enters any derivation; measured values appear only in the comparison columns.

## Requirements

- Python 3
- NumPy 2.0 or later — the scripts use `numpy.trapezoid`, which does not exist in NumPy 1.x
- SciPy 1.10 or later

```bash
pip install "numpy>=2.0" scipy
```

## Running

```bash
cd paper1
python physics_from_existence_pure.py
python physics_from_existence_verify.py
```

The two scripts separate derivation from verification:

| Script | What it does |
| --- | --- |
| `physics_from_existence_pure.py` | Derives the physical quantities from $V = -H$ alone. No experimental value appears anywhere in the code; the comments follow the paper section by section, so the derivation can be read independently. |
| `physics_from_existence_verify.py` | Runs the same derivation, then compares the results with experiment (CODATA 2022, PDG 2026, NuFIT 6.1). The experimental values are loaded in a separate function and used only for comparison — never as input. |

Both exit with status 0 on success. The verification script prints the paper version it corresponds to in its banner; check that it matches the [current release](../04-releases/README.md).

## Reading the output

The verification script works through the paper's chain in order:

1. **The equation and its bound states.** It solves the Schrödinger equation for $V = -H(\sigma(\phi))$ on a grid and prints the three bound-state energies together with the first continuum level, so the count $N = 3$ is visible directly. The WKB count and the quantity $\varepsilon$ derived from it follow, alongside the independently certified interval values for comparison.
2. **Derived quantities.** Gauge couplings, the Higgs ratio, the mass ratios, the mixing angles and the neutrino masses, each printed next to the measured value and the relative difference.
3. **Coverage table.** At the end, every row of the paper's Tables 1–4 is listed with one of four labels:

| Label | Meaning |
| --- | --- |
| evaluated | Computed by this script and compared with experiment |
| certified | Enclosed by interval arithmetic (the certificates are separate scripts, not part of this repository) |
| comparison-only | The measured value is shown; the derivation is in Paper II |
| not implemented (Paper II) | Derived in Paper II and not recomputed here |

No row is claimed as reproduced unless it is marked *evaluated* or *certified*. This is the honest boundary of what the repository verifies on its own.

## Precision

The grid solver works in floating point; the displayed digits of the paper come from the higher-precision and interval-certified evaluations described in the paper, and the script prints both where they differ (for example the WKB count, where the grid value and the certified value agree to the displayed digits). Where the paper quotes a value "at all orders", the script prints the closed-form evaluation and the certified interval it falls in.
