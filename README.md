# Physics from Existence

**Author:** Dr. Hidekazu Kondo, Kondo Research Institute (Tokyo, Japan)

Formalizing what it means "to exist" uniquely yields a single equation, *V = −H*.
This equation contains no free parameters, yet derives 26 fundamental constants of physics.
This suggests the existence of a new framework in which information theory describes matter and the universe.

**Current version:** v2.5.0 — *Physics from Existence I: The Equation*

**Web:** [Kondo Research Institute](https://www.kondo-lab.com/)

**Zenodo:** https://doi.org/10.5281/zenodo.19451887 (all versions)

## Comparison with Experiment

| Parameter | PFE | Experiment | Accuracy |
|-----------|-----|------------|----------|
| 1/α_em | 137.0368 | 137.036 | 0.0006% |
| sin θ_C | 0.2245 | 0.2243 | 0.10% |
| R_lepton | 1.5293 | 1.5294 | 0.006% |
| m_H/v | 0.5084 | 0.5082 | 0.03% |
| α_s | 0.1179 | 0.1179 | 0.01% |
| sin²θ_W | 3/13 | 0.23122 | 0.2% |
| R_up | 1.777 | 1.770 | 0.4% |
| R_down | 2.273 | 2.276 | 0.1% |
| sin²θ₁₂ | 4/13 | 0.3088 | 0.4% |
| sin²θ₁₃ | 0.0219 | 0.02249 | 2.6% |
| sin²θ₂₃ | 0.5475 | 0.561 | 2.4% |
| m₁ | 0.32 meV | — | prediction |
| θ_QCD | 0 | < 10⁻¹⁰ | exact |

*Current deviations are due to numerical precision and can be improved with higher-order calculations.*

## Version History

| Version | Zenodo | Date | Title | Zenodo DOI |
|---------|--------|------|-------|------------|
| v2.5.0 | v12 | 2026-07-23 | Physics from Existence I: The Equation | [10.5281/zenodo.21511757](https://zenodo.org/records/21511757) |
| v2.4.0 | v11 | 2026-07-02 | Physics from Existence I: The Equation | [10.5281/zenodo.21123893](https://zenodo.org/records/21123893) |
| v2.3.0 | v10 | 2026-06-12 | Physics from Existence I: The Equation | [10.5281/zenodo.20666890](https://zenodo.org/records/20666890) |
| v2.2.1 | v9 | 2026-06-11 | Physics from Existence I: The Equation | [10.5281/zenodo.20634925](https://zenodo.org/records/20634925) |
| v2.2.0 | v8 | 2026-06-02 | Physics from Existence I: The Equation | [10.5281/zenodo.20488132](https://zenodo.org/records/20488132) |
| v2.1.0 | v7 | 2026-05-31 | Physics from Existence I: The Equation | [10.5281/zenodo.20471724](https://zenodo.org/records/20471724) |
| v2.0.0 | v6 | 2026-05-31 | Physics from Existence I: The Equation | [10.5281/zenodo.20470331](https://zenodo.org/records/20470331) |
| v1.2.0 | v5 | 2026-04-21 | Physics from Existence I: The Equation | [10.5281/zenodo.19678855](https://doi.org/10.5281/zenodo.19678855) |
| v1.1.0 | v4 | 2026-04-13 | Physics from Existence I: The Equation | [10.5281/zenodo.19551657](https://doi.org/10.5281/zenodo.19551657) |
| v1.0.2 | v3 | 2026-04-10 | Physics from Existence | [10.5281/zenodo.19488819](https://doi.org/10.5281/zenodo.19488819) |
| v1.0.1 | v2 | 2026-04-09 | Physics from Existence | [10.5281/zenodo.19478719](https://doi.org/10.5281/zenodo.19478719) |
| v1.0.0 | v1 | 2026-04-07 | Physics from Existence | [10.5281/zenodo.19451888](https://doi.org/10.5281/zenodo.19451888) |

## Repository Structure

```
paper1/
  physics_from_existence.tex         LaTeX source (English)
  physics_from_existence.pdf         Compiled PDF (English)
  physics_from_existence_ja.tex      LaTeX source (Japanese)
  physics_from_existence_ja.pdf      Compiled PDF (Japanese)
  physics_from_existence_pure.py     Pure derivation (no experimental input)
  physics_from_existence_verify.py   Verification against experiment
```

## Scripts

Two scripts are provided to separate derivation from verification:

**`pure.py`** — Derives 20+ physical quantities from `V = -H` alone. No experimental values appear anywhere in the code. Comments follow the paper section by section, making the derivation independently readable.

```bash
cd paper1 && python physics_from_existence_pure.py
```

**`verify.py`** — Runs the same derivation, then compares results against experimental data (CODATA, PDG, NuFIT). The experimental values are loaded in a separate function (`load_experiment`) and used only for comparison — never as input to the derivation.

```bash
cd paper1 && python physics_from_existence_verify.py
```

Both require Python 3, NumPy, and SciPy.

## Formats

All content is stored in plain-text formats (LaTeX, Markdown, Python) to
ensure readability without specialized software, indefinitely.

## License

[CC BY 4.0](LICENSE) — free to share and adapt with attribution.
