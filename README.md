# Physics from Existence

**Author:** Dr. Hidekazu Kondo, Kondo Research Institute (Tokyo, Japan)

Formalizing what it means "to exist" uniquely yields a single equation, *V = −H*.
This equation contains no free parameters, yet the 26 physical constants of the Standard Model are reproduced from it with no adjustable parameters.
This suggests the existence of a new framework in which information theory describes matter and the universe.

**Current version:** v4.1.1 — *Physics from Existence I: The Equation* (2026-09-01)

**Zenodo:** https://zenodo.org/records/22217837 · all versions: https://doi.org/10.5281/zenodo.19451887

**Web:** [Kondo Research Institute](https://www.kondo-lab.com/)

## Documentation

The [`docs/`](docs/README.md) folder explains the repository:

| | |
| --- | --- |
| [Papers](docs/01-papers/README.md) | Abstract, files, how to cite |
| [Verification](docs/02-verification/README.md) | Running the scripts and reading their output |
| [Results](docs/03-results/README.md) | Comparison with the Standard Model constants, structural results, testable predictions |
| [Releases](docs/04-releases/README.md) | Version history with Zenodo DOIs, versioning policy, licences |

日本語版は [`docs/i18n/ja/`](docs/i18n/ja/README.md) にあります。

## Repository Structure

```
paper1/
  physics_from_existence.tex         LaTeX source (English)
  physics_from_existence.pdf         Compiled PDF (English)
  physics_from_existence_ja.tex      LaTeX source (Japanese)
  physics_from_existence_ja.pdf      Compiled PDF (Japanese)
  physics_from_existence_pure.py     Pure derivation (no experimental input)
  physics_from_existence_verify.py   Verification against experiment
docs/                                Documentation (English, with Japanese under i18n/ja/)
```

File names are fixed across versions; each release replaces their contents and is tagged `vX.Y.Z`.

## Quick verification

```bash
pip install "numpy>=2.0" scipy
cd paper1 && python physics_from_existence_verify.py
```

`pure.py` derives 20+ physical quantities from *V = −H* alone, with no experimental value anywhere in the code. `verify.py` runs the same derivation and compares the results with CODATA, PDG and NuFIT values, which are loaded separately and used only for comparison. See [Verification](docs/02-verification/README.md).

## Formats

All content is stored in plain-text formats (LaTeX, Markdown, Python) to
ensure readability without specialized software, indefinitely.

## License

- **Paper, text, and figures** (`.tex`, `.pdf`, figures) — [CC BY 4.0](LICENSE): free to share and adapt with attribution.
- **Code** (`.py`) — [MIT](LICENSE-CODE): free to use, modify, and redistribute, including within other projects.
