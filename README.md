# Physics from Existence

**Author:** Dr. Hidekazu Kondo (Tokyo, Japan)

Formalizing what it means "to exist" uniquely yields a single equation, *V = −H*.
This equation contains no free parameters, yet derives 26 fundamental constants of physics.
This suggests the existence of a new framework in which information theory describes matter and the universe.

**Current version:** v2.1.0 — *Physics from Existence I: The Equation*

**Zenodo:** https://doi.org/10.5281/zenodo.19451887 (all versions)

**Latest:** https://zenodo.org/records/20471724

## Version History

| Version | Date | Title | Zenodo DOI |
|---------|------|-------|------------|
| v2.1.0 | 2026-05-31 | Physics from Existence I: The Equation | [10.5281/zenodo.20471724](https://zenodo.org/records/20471724) |
| v1.0.2 | 2026-04-10 | Physics from Existence | [10.5281/zenodo.19488819](https://doi.org/10.5281/zenodo.19488819) |
| v1.0.1 | 2026-04-09 | Physics from Existence | [10.5281/zenodo.19478719](https://doi.org/10.5281/zenodo.19478719) |
| v1.0.0 | 2026-04-07 | Physics from Existence | [10.5281/zenodo.19451888](https://doi.org/10.5281/zenodo.19451888) |

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
