# Physics from Existence

**Author:** Dr. Hidekazu Kondo (Tokyo, Japan)

Formalizing what it means "to exist" uniquely yields a single equation, *V = −H*.
This equation contains no free parameters, yet derives 26 fundamental constants of physics.
This suggests the existence of a new framework in which information theory describes matter and the universe.

**Current version:** v2.2.0 — *Physics from Existence I: The Equation*

**Zenodo:** https://doi.org/10.5281/zenodo.19451887 (all versions)

## Comparison with Experiment

| Parameter | PFE | Experiment | Accuracy |
|-----------|-----|------------|----------|
| 1/α_em | 137.0367 | 137.036 | 0.0006% |
| sin θ_C | 0.2245 | 0.2244 | 0.05% |
| R_lepton | 1.5293 | 1.5294 | 0.006% |
| m_H/v | 0.5084 | 0.5087 | 0.07% |
| α_s | 0.1179 | 0.1180 | 0.08% |
| sin²θ_W | 3/13 | 0.23122 | 0.2% |
| R_up | 1.777 | 1.772 | 0.3% |
| R_down | 2.273 | 2.269 | 0.2% |
| sin²θ₁₂ | 4/13 | 0.307 | 0.2% |
| sin²θ₁₃ | 0.0219 | 0.02195 | 0.3% |
| sin²θ₂₃ | 0.5475 | 0.561 | 2.4% |
| m₁ | 0.31 meV | — | prediction |
| θ_QCD | 0 | < 10⁻¹⁰ | exact |

## Version History

| Version | Zenodo | Date | Title | Zenodo DOI |
|---------|--------|------|-------|------------|
| v2.2.0 | v8 | 2026-06-02 | Physics from Existence I: The Equation | [10.5281/zenodo.20488132](https://zenodo.org/records/20488132) |
| v2.1.0 | v7 | 2026-05-31 | Physics from Existence I: The Equation | [10.5281/zenodo.20471724](https://zenodo.org/records/20471724) |
| v1.0.2 | v3 | 2026-04-10 | Physics from Existence | [10.5281/zenodo.19488819](https://doi.org/10.5281/zenodo.19488819) |
| v1.0.1 | v2 | 2026-04-09 | Physics from Existence | [10.5281/zenodo.19478719](https://doi.org/10.5281/zenodo.19478719) |
| v1.0.0 | v1 | 2026-04-07 | Physics from Existence | [10.5281/zenodo.19451888](https://doi.org/10.5281/zenodo.19451888) |

## Repository Structure

```
paper1/
  physics_from_existence.tex         LaTeX source
  physics_from_existence.pdf         Compiled PDF
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
