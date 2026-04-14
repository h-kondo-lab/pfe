# Physics from Existence

**Author:** Hidekazu Kondo (Tokyo, Japan)

The canonical free energy of a single fermionic mode defines a parameter-free
potential whose spectral properties reproduce the Standard Model's dimensionless
constants — including 1/alpha = 137.037, three generations, and the gauge
structure — with zero free parameters.

**Zenodo:** https://doi.org/10.5281/zenodo.19451887 (all versions)

## Version History

| Version | Date | Zenodo DOI |
|---------|------|------------|
| v1.0.0 | 2026-04-07 | [10.5281/zenodo.19451888](https://doi.org/10.5281/zenodo.19451888) |
| v1.0.1 | 2026-04-09 | [10.5281/zenodo.19478719](https://doi.org/10.5281/zenodo.19478719) |
| v1.0.2 | 2026-04-10 | [10.5281/zenodo.19488819](https://doi.org/10.5281/zenodo.19488819) |

## Repository Structure

```
paper1/
  physics_from_existence.tex      LaTeX source
  physics_from_existence.pdf      Compiled PDF
  physics_from_existence_verify.py  Numerical verification script
```

## Verification

```bash
cd paper1 && python physics_from_existence_verify.py
```

All numerical predictions are independently computed.
No dependencies beyond Python 3 and NumPy.

## Formats

All content is stored in plain-text formats (LaTeX, Markdown, Python) to
ensure readability without specialized software, indefinitely.

## License

[CC BY 4.0](LICENSE) — free to share and adapt with attribution.
